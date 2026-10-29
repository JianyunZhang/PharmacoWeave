"""Failure and omission rules for the locked IV bolus NCA calculator."""

import math

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import ObservationStatus, Route
from pharmacoweave.errors import UnsupportedDesign
from pharmacoweave.pk_nca.calculator import calculate_profile
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_nca.runner import run_nca
from pharmacoweave.pk_profiles import ConcentrationProfile, TimedConcentration


def _config(**overrides: object) -> NcaConfig:
    payload: dict[str, object] = {
        "analysis": "nca",
        "design_route": "iv_bolus",
        "dosing": "single_dose",
        "sampling": "rich_serial",
        "auc_method": "linear_up_log_down",
        "lambda_z_min_points": 3,
        "lambda_z_min_adjusted_r2": 0.90,
        "auc_inf_max_extrapolated_fraction": 0.20,
        "blq_rule": "keep_as_censored",
    }
    payload.update(overrides)
    return NcaConfig.model_validate(payload)


def _profile(
    points: list[tuple[float, ObservationStatus, float | None]],
    *,
    dose: float = 1.0,
    dose_unit: str = "mg",
    unit: str = "mg/L",
) -> ConcentrationProfile:
    timed = tuple(
        TimedConcentration(
            time_h=time_h,
            status=status,
            value=value,
            value_unit=unit if status is ObservationStatus.QUANTIFIED else None,
            observation_id=f"OBS-{index}",
            is_predose=time_h < 0,
        )
        for index, (time_h, status, value) in enumerate(points)
    )
    return ConcentrationProfile(
        subject_id="SUBJ-1",
        analyte_id="ANALYTE-1",
        dose_event_id="DOSE-1",
        route=Route.IV_BOLUS,
        dose_value=dose,
        dose_unit=dose_unit,
        concentration_unit=unit,
        lloq=None,
        points=timed,
    )


def test_linear_up_and_log_down_are_both_used() -> None:
    record = calculate_profile(
        _profile(
            [
                (0.0, ObservationStatus.QUANTIFIED, 1.0),
                (1.0, ObservationStatus.QUANTIFIED, 3.0),
                (2.0, ObservationStatus.QUANTIFIED, 1.0),
            ]
        ),
        _config(),
    )
    assert record.auc_last == pytest.approx(2.0 + 2.0 / math.log(3.0))
    assert record.cmax == pytest.approx(3.0)
    assert record.tmax == pytest.approx(1.0)
    assert record.terminal_phase_qualified is False
    assert record.lambda_z is None
    assert record.auc_inf is None


def test_exact_monoexponential_tie_uses_the_longer_terminal_regression() -> None:
    points = [(float(time), ObservationStatus.QUANTIFIED, math.exp(-0.2 * time)) for time in range(1, 7)]
    record = calculate_profile(_profile(points, dose=1.0), _config())
    assert record.terminal_phase_n_points == 5
    assert record.lambda_z == pytest.approx(0.2, abs=1e-9)
    assert record.terminal_phase_qualified is True


def test_too_few_points_omit_the_terminal_phase(repo_root) -> None:
    study = _filter_minutes(build_study(), {-5, 5, 15})
    records = run_nca(study, repo_root / "configs" / "nca_single_dose_iv.yaml")
    assert records.nca_parameters
    for record in records.nca_parameters:
        assert record.cmax is not None
        assert record.auc_last is not None
        assert record.terminal_phase_qualified is False
        assert record.lambda_z is None
        assert record.half_life is None
        assert record.auc_inf is None


def test_low_adjusted_r2_omits_lambda_and_auc_infinity() -> None:
    record = calculate_profile(
        _profile(
            [
                (0.1, ObservationStatus.QUANTIFIED, 8.0),
                (1.0, ObservationStatus.QUANTIFIED, 5.0),
                (2.0, ObservationStatus.QUANTIFIED, 4.5),
                (3.0, ObservationStatus.QUANTIFIED, 1.2),
                (4.0, ObservationStatus.QUANTIFIED, 3.8),
                (5.0, ObservationStatus.QUANTIFIED, 0.4),
            ]
        ),
        _config(),
    )
    assert record.terminal_phase_qualified is False
    assert record.lambda_z is None
    assert record.auc_inf is None
    assert record.auc_last is not None
    assert any("not qualified" in note for note in record.warnings)


def test_large_extrapolation_omits_auc_infinity_and_keeps_lambda() -> None:
    points = []
    for minute in (5, 15, 30, 60):
        time_h = minute / 60
        points.append((time_h, ObservationStatus.QUANTIFIED, round(math.exp(-0.2 * time_h), 6)))
    record = calculate_profile(_profile(points, dose=0.25), _config())
    assert record.terminal_phase_qualified is True
    assert record.lambda_z == pytest.approx(0.2, abs=1e-4)
    assert record.auc_inf is None
    assert record.clearance is None
    assert record.volume is None
    assert any("above the reporting limit" in note for note in record.warnings)


def test_internal_blq_is_not_imputed_or_bridged() -> None:
    record = calculate_profile(
        _profile(
            [
                (0.0, ObservationStatus.QUANTIFIED, 1.0),
                (1.0, ObservationStatus.BLQ, None),
                (2.0, ObservationStatus.QUANTIFIED, 0.5),
            ]
        ),
        _config(),
    )
    assert record.auc_last is None
    assert record.cmax == pytest.approx(1.0)
    assert record.excluded_reason is not None
    assert "censored" in record.excluded_reason
    excluded = calculate_profile(
        _profile(
            [
                (0.0, ObservationStatus.QUANTIFIED, 1.0),
                (1.0, ObservationStatus.BLQ, None),
                (2.0, ObservationStatus.QUANTIFIED, 0.5),
            ]
        ),
        _config(blq_rule="exclude"),
    )
    assert excluded.auc_last == pytest.approx(1.0 / math.log(2.0))
    assert any("bridges" in note for note in excluded.warnings)


def test_duplicate_times_and_incompatible_dose_units_omit_parameters() -> None:
    duplicate = calculate_profile(
        _profile(
            [
                (1.0, ObservationStatus.QUANTIFIED, 1.0),
                (1.0, ObservationStatus.QUANTIFIED, 0.9),
                (2.0, ObservationStatus.QUANTIFIED, 0.5),
            ]
        ),
        _config(),
    )
    assert duplicate.auc_last is None
    assert duplicate.excluded_reason is not None
    assert "Duplicate" in duplicate.excluded_reason
    mismatched = calculate_profile(
        _profile(
            [
                (float(time), ObservationStatus.QUANTIFIED, math.exp(-0.3 * time))
                for time in (0, 1, 2, 4, 6, 8, 10)
            ],
            dose=1.0,
            dose_unit="nmol",
            unit="mg/L",
        ),
        _config(),
    )
    assert mismatched.auc_last is not None
    assert mismatched.clearance is None
    assert any("Clearance" in note for note in mismatched.warnings)


def test_oral_multiple_dose_and_sparse_designs_are_rejected(repo_root) -> None:
    study = build_study()
    config = repo_root / "configs" / "nca_single_dose_iv.yaml"
    oral = study.model_copy(
        update={"dose_events": [dose.model_copy(update={"route": Route.ORAL}) for dose in study.dose_events]}
    )
    with pytest.raises(UnsupportedDesign, match="iv_bolus"):
        run_nca(oral, config)
    extra = study.dose_events[0].model_copy(update={"dose_event_id": "DOSE-EXTRA"})
    repeated = study.model_copy(update={"dose_events": [*study.dose_events, extra]})
    with pytest.raises(UnsupportedDesign, match="single_dose"):
        run_nca(repeated, config)
    sparse = _filter_minutes(study, {-5, 5})
    with pytest.raises(UnsupportedDesign, match="Sparse"):
        run_nca(sparse, config)


def _filter_minutes(study, minutes: set[float]):
    samples = [sample for sample in study.biosamples if sample.time_since_dose in minutes]
    identifiers = {sample.sample_id for sample in samples}
    observations = [item for item in study.observations if item.sample_id in identifiers]
    return study.model_copy(update={"biosamples": samples, "observations": observations})
