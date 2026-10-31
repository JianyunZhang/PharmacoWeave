"""Compare quiet single-dose constant-rate IV infusion NCA with the frozen 50-digit reference."""

import json
from pathlib import Path

import pytest
from examples.iv_infusion_case import build_iv_infusion_study

from pharmacoweave.core.entities import Route
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint, iv_infusion_profiles
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_nca.calculate import nca_iv_infusion
from pharmacoweave.pk_nca.config import NcaConfig

_REFERENCE = Path(__file__).resolve().parent / "iv_infusion_quiet_nca.json"
_DOCUMENTED_RELATIVE = {
    "auc_last": 1e-9,
    "auc_inf": 1e-9,
    "auc_inf_extrap_fraction": 1e-9,
    "lambda_z": 1e-9,
    "half_life": 1e-9,
    "clearance": 1e-9,
    "volume": 1e-9,
    "terminal_phase_adjusted_r2": 1e-9,
    "cmax": 1e-12,
}
_DOCUMENTED_ABSOLUTE = {"tmax_h": 1e-12}


def _load() -> dict:
    return json.loads(_REFERENCE.read_text(encoding="utf-8"))


def test_iv_infusion_acceptance_bounds_match_the_validation_document() -> None:
    reference = _load()
    assert reference["acceptance"]["relative"] == _DOCUMENTED_RELATIVE
    assert reference["acceptance"]["absolute"] == _DOCUMENTED_ABSOLUTE
    assert "not produced by importing" in reference["computed_with"]
    omitted = reference["omitted_observations"]
    assert any(item["time_min"] == "1440" and item["status"] == "BLQ" for item in omitted)
    assert all(point["concentration"] != "0" for point in reference["quantified_concentrations_mg_l"])


def test_iv_infusion_profiles_use_the_frozen_quantified_points() -> None:
    reference = _load()
    expected = reference["quantified_concentrations_mg_l"]
    for profile in iv_infusion_profiles(build_iv_infusion_study()):
        assert profile.route is Route.IV_INFUSION
        assert profile.infusion_duration_h == pytest.approx(1.0, abs=1e-12)
        assert profile.n_blq_omitted == 1
        assert profile.n_predose_omitted == 1
        assert len(profile.points) == len(expected)
        for point, row in zip(profile.points, expected, strict=True):
            assert point.time_h == pytest.approx(float(row["time_h"]), abs=1e-12)
            assert point.concentration_mg_l == pytest.approx(float(row["concentration"]), rel=1e-12)
            assert point.concentration_mg_l > 0


def test_iv_infusion_nca_matches_the_frozen_reference() -> None:
    reference = _load()
    shared = reference["shared_parameters"]
    relative = reference["acceptance"]["relative"]
    absolute = reference["acceptance"]["absolute"]
    by_subject = {
        row.subject_id: row
        for row in build_iv_infusion_study().run_nca("configs/nca_single_dose_iv_infusion.yaml")
    }
    subjects = {item["subject_id"]: item for item in reference["subjects"]}
    assert set(by_subject) == set(subjects)
    for subject_id, row in by_subject.items():
        frozen = subjects[subject_id]
        assert row.terminal_phase_qualified
        assert row.terminal_phase_n_points == shared["terminal_phase_n_points"]
        assert row.auc_last == pytest.approx(float(shared["auc_last_mg_h_l"]), rel=relative["auc_last"])
        assert row.auc_inf == pytest.approx(float(shared["auc_inf_mg_h_l"]), rel=relative["auc_inf"])
        assert row.auc_inf_extrap_fraction == pytest.approx(
            float(shared["auc_inf_extrap_fraction"]),
            rel=relative["auc_inf_extrap_fraction"],
        )
        assert row.lambda_z == pytest.approx(float(shared["lambda_z_per_h"]), rel=relative["lambda_z"])
        assert row.half_life == pytest.approx(float(shared["half_life_h"]), rel=relative["half_life"])
        assert row.cmax == pytest.approx(float(shared["cmax_mg_l"]), rel=relative["cmax"])
        assert row.tmax == pytest.approx(float(shared["tmax_h"]), abs=absolute["tmax_h"])
        assert row.terminal_phase_adjusted_r2 == pytest.approx(
            float(shared["terminal_phase_adjusted_r2"]),
            rel=relative["terminal_phase_adjusted_r2"],
        )
        assert row.clearance == pytest.approx(float(frozen["clearance_l_h"]), rel=relative["clearance"])
        assert row.volume == pytest.approx(float(frozen["volume_l"]), rel=relative["volume"])
        assert row.clearance_semantic is not None and row.clearance_semantic.value == "CL"
        assert row.volume_semantic is not None and row.volume_semantic.value == "Vz"
        assert not any("C0 was back-extrapolated" in warning for warning in row.warnings)
        assert any("not replaced with a number" in warning for warning in row.warnings)


def test_iv_infusion_nca_requires_positive_infusion_duration() -> None:
    cfg = NcaConfig.model_validate(
        {
            "analysis": "nca",
            "design_route": "iv_infusion",
            "dosing": "single_dose",
            "sampling": "rich_serial",
            "auc_method": "linear_up_log_down",
            "lambda_z_min_points": 3,
            "lambda_z_min_adjusted_r2": 0.90,
            "auc_inf_max_extrapolated_fraction": 0.20,
            "blq_rule": "keep_as_censored",
        }
    )
    missing_duration = ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="PW-0001-PARENT",
        dose_event_id="DOSE-1",
        route=Route.IV_INFUSION,
        dose_mg=0.5,
        points=(
            QuantifiedPoint("OBS-1", 1.0, 1.5),
            QuantifiedPoint("OBS-2", 2.0, 1.1),
            QuantifiedPoint("OBS-3", 4.0, 0.6),
        ),
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
        infusion_duration_h=None,
    )
    with pytest.raises(AnalysisError, match="infusion_duration"):
        nca_iv_infusion(missing_duration, cfg)
