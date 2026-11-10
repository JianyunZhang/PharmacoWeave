"""Phase 4 slices from the written plan. NLME, dynamic PBPK, and QSP networks stay refused."""

import ast
import json
import math
import subprocess
from pathlib import Path

import pytest
from pydantic import ValidationError

from pharmacoweave.errors import AnalysisError, PlannedCapabilityError
from pharmacoweave.pbpk import run_in_process as run_pbpk
from pharmacoweave.pbpk.physiology_table import PhysiologyTable, TissuePhysiology, validate_physiology_table
from pharmacoweave.pbpk.steady_state import steady_state_one_tissue
from pharmacoweave.population_pk import run_in_process as run_population
from pharmacoweave.population_pk.results import load_population_result_table
from pharmacoweave.population_pk.two_stage import IndividualParameterEstimate, summarize_log_parameters
from pharmacoweave.qsp import run_in_process as run_qsp
from pharmacoweave.qsp.turnover import simulate_linear_turnover

_REPO = Path(__file__).resolve().parents[2]
_TWO_STAGE = json.loads((_REPO / "tests" / "regression" / "two_stage_log_parameters.json").read_text(encoding="utf-8"))
_UNITS = json.loads((_REPO / "tests" / "regression" / "two_stage_canonical_units.json").read_text(encoding="utf-8"))
_FLOW = json.loads((_REPO / "tests" / "regression" / "physiology_flow_tolerance.json").read_text(encoding="utf-8"))
_PBPK = json.loads((_REPO / "tests" / "regression" / "steady_state_one_tissue.json").read_text(encoding="utf-8"))
_TURNOVER = json.loads((_REPO / "tests" / "regression" / "linear_turnover.json").read_text(encoding="utf-8"))


def test_two_stage_matches_the_hand_formula_and_rejects_bad_inputs() -> None:
    estimates = [
        _estimate(f"SUBJ-0{index}", "CL", value)
        for index, value in enumerate(_TWO_STAGE["estimates_l_h"], start=1)
    ]
    summary = summarize_log_parameters(estimates)[0]
    hand_mean = math.exp(math.fsum(math.log(value) for value in _TWO_STAGE["estimates_l_h"]) / 2)
    hand_variance = math.fsum(
        (math.log(value) - math.log(_TWO_STAGE["geometric_mean_l_h"])) ** 2 for value in _TWO_STAGE["estimates_l_h"]
    )
    assert summary.geometric_mean == pytest.approx(_TWO_STAGE["geometric_mean_l_h"], rel=1e-12, abs=1e-12)
    assert summary.geometric_mean == pytest.approx(hand_mean, rel=1e-12, abs=1e-12)
    assert summary.log_variance == pytest.approx(hand_variance, rel=1e-12, abs=1e-12)
    assert summary.n_subjects == 2
    assert any("not an NLME posterior" in item for item in summary.limitations)

    with pytest.raises(AnalysisError, match="fewer than two"):
        summarize_log_parameters([_estimate("SUBJ-01", "CL", 0.05)])
    with pytest.raises(AnalysisError, match="cannot be summarized"):
        summarize_log_parameters([_estimate("SUBJ-01", "CL", 0.05), _estimate("SUBJ-02", "CL/F", 0.20)])
    with pytest.raises(AnalysisError, match="positive"):
        summarize_log_parameters([_estimate("SUBJ-01", "CL", 0.05), _estimate("SUBJ-02", "CL", 0.0)])
    with pytest.raises(AnalysisError, match="repeats"):
        summarize_log_parameters([_estimate("SUBJ-01", "CL", 0.05), _estimate("SUBJ-01", "CL", 0.20)])
    wrong_unit = IndividualParameterEstimate(subject_id="SUBJ-01", parameter_name="CL", estimate=0.05, unit="mg")
    with pytest.raises(AnalysisError, match="dimension"):
        summarize_log_parameters([wrong_unit, _estimate("SUBJ-02", "CL", 0.20)])


def test_physiology_check_has_no_concentration_and_rejects_unbalanced_flow() -> None:
    checked = validate_physiology_table(_physiology())
    assert checked.passed is True
    dumped = checked.model_dump()
    assert "concentration" not in dumped
    assert "auc" not in dumped
    assert "clearance" not in dumped
    with pytest.raises(AnalysisError, match="does not match cardiac output"):
        validate_physiology_table(_physiology(muscle_flow=1.5))
    duplicate = _physiology()
    duplicate.tissues[1] = duplicate.tissues[0].model_copy(update={"tissue_name": "liver"})
    with pytest.raises(AnalysisError, match="once"):
        validate_physiology_table(duplicate)
    with pytest.raises(ValidationError):
        TissuePhysiology(
            tissue_name="liver",
            volume=0,
            volume_unit="L",
            blood_flow=1,
            blood_flow_unit="L/h",
            partition_coefficient=1,
            partition_unit="1",
            source="caller note",
        )
    with pytest.raises(ValidationError):
        TissuePhysiology(
            tissue_name="liver",
            volume=1,
            volume_unit="L",
            blood_flow=1,
            blood_flow_unit="L/h",
            partition_coefficient=1,
            partition_unit="1",
            source="x",
        )


def test_turnover_matches_the_closed_form() -> None:
    result = simulate_linear_turnover(
        kin=_TURNOVER["kin"],
        kout=_TURNOVER["kout"],
        r0=_TURNOVER["r0"],
        times_h=list(_TURNOVER["times_h"]),
        context_of_use="Research check of one linear turnover equation.",
    )
    assert result.evidence_grade == "predicted"
    assert result.steady_state == pytest.approx(_TURNOVER["steady_state"], rel=1e-12, abs=1e-12)
    assert result.curve[0].response == pytest.approx(_TURNOVER["response_at_0"], abs=1e-12)
    expected_one = (_TURNOVER["kin"] / _TURNOVER["kout"]) * (1.0 - math.exp(-_TURNOVER["kout"] * 1.0))
    assert result.curve[1].response == pytest.approx(expected_one, rel=1e-12, abs=1e-12)
    assert result.constitutes_dosing_recommendation is False
    with pytest.raises(AnalysisError, match="positive"):
        simulate_linear_turnover(kin=0, kout=1, r0=0, times_h=[0, 1], context_of_use="research")
    with pytest.raises(AnalysisError, match="does not estimate"):
        simulate_linear_turnover(
            kin=2, kout=1, r0=0, times_h=[0, 1], context_of_use="research", estimate_parameters=True
        )
    with pytest.raises(AnalysisError, match="strictly increase"):
        simulate_linear_turnover(kin=2, kout=1, r0=0, times_h=[1, 1], context_of_use="research")


def test_one_tissue_steady_state_matches_the_hand_arithmetic() -> None:
    result = steady_state_one_tissue(
        _physiology(),
        input_rate=_PBPK["input_rate_mg_h"],
        input_rate_unit="mg/h",
        hepatic_intrinsic_clearance=_PBPK["hepatic_intrinsic_clearance_l_h"],
        clearance_unit="L/h",
        clearance_name="CL",
    )
    assert result.evidence_grade == "predicted"
    assert result.clearance_name == "CL"
    assert result.hepatic_clearance_l_h == pytest.approx(_PBPK["hepatic_clearance_l_h"], rel=1e-12, abs=1e-12)
    assert result.arterial_concentration_mg_l == pytest.approx(_PBPK["arterial_concentration_mg_l"], rel=1e-12)
    assert result.liver_concentration_mg_l == pytest.approx(_PBPK["liver_concentration_mg_l"], rel=1e-12)
    assert result.tissue_concentration_mg_l == pytest.approx(_PBPK["tissue_concentration_mg_l"], rel=1e-12)
    assert result.constitutes_dosing_recommendation is False
    with pytest.raises(AnalysisError, match="time course"):
        steady_state_one_tissue(
            _physiology(),
            input_rate=1,
            input_rate_unit="mg/h",
            hepatic_intrinsic_clearance=1,
            clearance_unit="L/h",
            clearance_name="CL/F",
            time_course=True,
        )


def test_extra_tissue_and_population_table_rules(monkeypatch: pytest.MonkeyPatch) -> None:
    kidney = TissuePhysiology(
        tissue_name="kidney",
        volume=0.2,
        volume_unit="L",
        blood_flow=0.5,
        blood_flow_unit="L/h",
        partition_coefficient=1,
        partition_unit="1",
        source="caller kidney note",
    )
    table = _physiology(muscle_flow=0.5)
    table.tissues.append(kidney)
    with pytest.raises(AnalysisError, match="exactly one other tissue"):
        steady_state_one_tissue(
            table,
            input_rate=1,
            input_rate_unit="mg/h",
            hepatic_intrinsic_clearance=1,
            clearance_unit="L/h",
            clearance_name="CL",
        )

    loaded = load_population_result_table(
        [
            {
                "subject_id": "SUBJ-01",
                "parameter_name": "CL/F",
                "estimate": 0.2,
                "unit": "L/h",
                "estimation_method": "caller",
                "evidence_grade": "estimated",
            }
        ]
    )
    assert loaded.rows[0].evidence_grade == "estimated"
    assert loaded.rows[0].parameter_name == "CL/F"

    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("reading a population table started a subprocess")

    monkeypatch.setattr(subprocess, "run", fail_if_called)
    with pytest.raises(AnalysisError, match="missing"):
        load_population_result_table([{"parameter_name": "CL", "estimate": 1}])
    with pytest.raises(AnalysisError, match="relabels"):
        load_population_result_table(
            [
                {
                    "subject_id": "SUBJ-01",
                    "parameter_name": "CL",
                    "source_parameter_name": "CL/F",
                    "estimate": 0.2,
                    "unit": "L/h",
                    "estimation_method": "caller",
                    "evidence_grade": "estimated",
                }
            ]
        )
    with pytest.raises(AnalysisError, match="evidence grade"):
        load_population_result_table(
            [
                {
                    "subject_id": "SUBJ-01",
                    "parameter_name": "CL",
                    "estimate": 0.2,
                    "unit": "L/h",
                    "estimation_method": "caller",
                    "evidence_grade": "observed",
                }
            ]
        )
    with pytest.raises(AnalysisError, match="dosing recommendation"):
        load_population_result_table(
            [
                {
                    "subject_id": "SUBJ-01",
                    "parameter_name": "CL",
                    "estimate": 0.2,
                    "unit": "L/h",
                    "estimation_method": "caller",
                    "evidence_grade": "predicted",
                    "dosing_recommendation": 10,
                }
            ]
        )


def test_two_stage_drops_every_mean_when_one_group_fails_and_ignores_standard_error() -> None:
    quiet = summarize_log_parameters(
        [_estimate("SUBJ-01", "CL", 0.05), _estimate("SUBJ-02", "CL", 0.20)]
    )[0]
    weighted = summarize_log_parameters(
        [
            _estimate("SUBJ-01", "CL", 0.05, standard_error=1e-6),
            _estimate("SUBJ-02", "CL", 0.20, standard_error=50.0),
        ]
    )[0]
    assert weighted.geometric_mean == quiet.geometric_mean
    assert weighted.log_variance == quiet.log_variance
    assert "standard_error" not in weighted.model_dump()
    with pytest.raises(AnalysisError, match="cannot be summarized"):
        summarize_log_parameters([_estimate("SUBJ-01", "V", 1.0), _estimate("SUBJ-02", "V/F", 2.0)])
    with pytest.raises(AnalysisError, match="fewer than two"):
        summarize_log_parameters(
            [
                _estimate("SUBJ-01", "CL", 0.05),
                _estimate("SUBJ-02", "CL", 0.20),
                _estimate("SUBJ-01", "V", 1.0),
            ]
        )
    with pytest.raises(ValidationError):
        IndividualParameterEstimate.model_validate(
            {"subject_id": "SUBJ-01", "parameter_name": "CL", "unit": "L/h"}
        )


def test_two_stage_converts_clearance_and_volume_to_canonical_units() -> None:
    summaries = {
        item.parameter_name: item
        for item in summarize_log_parameters(
            [
                _estimate("SUBJ-01", "CL", _UNITS["clearance_estimates"][0], unit=_UNITS["clearance_unit"]),
                _estimate("SUBJ-02", "CL", _UNITS["clearance_estimates"][1], unit=_UNITS["clearance_unit"]),
                _estimate("SUBJ-01", "V", _UNITS["volume_estimates"][0], unit=_UNITS["volume_unit"]),
                _estimate("SUBJ-02", "V", _UNITS["volume_estimates"][1], unit=_UNITS["volume_unit"]),
            ]
        )
    }
    hand_variance = 2.0 * (math.log(2.0) ** 2)
    clearance = summaries["CL"]
    volume = summaries["V"]
    assert clearance.unit == "L/h"
    assert volume.unit == "L"
    assert clearance.geometric_mean == pytest.approx(_UNITS["clearance_geometric_mean_l_h"], rel=1e-12, abs=1e-12)
    assert volume.geometric_mean == pytest.approx(_UNITS["volume_geometric_mean_l"], rel=1e-12, abs=1e-12)
    assert clearance.log_variance == pytest.approx(hand_variance, rel=1e-12, abs=1e-12)
    assert volume.log_variance == pytest.approx(hand_variance, rel=1e-12, abs=1e-12)
    hand_clearance = math.exp(
        math.fsum(math.log(value) for value in _UNITS["clearance_canonical_l_h"]) / 2
    )
    hand_volume = math.exp(math.fsum(math.log(value) for value in _UNITS["volume_canonical_l"]) / 2)
    assert clearance.geometric_mean == pytest.approx(hand_clearance, rel=1e-12, abs=1e-12)
    assert volume.geometric_mean == pytest.approx(hand_volume, rel=1e-12, abs=1e-12)


def test_physiology_flow_tolerance_uses_the_frozen_band() -> None:
    relative = _FLOW["relative_tolerance"]
    absolute = _FLOW["absolute_tolerance"]
    cardiac = 2.0
    relative_limit = max(relative * cardiac, absolute)
    assert validate_physiology_table(_physiology(muscle_flow=1.0 + 0.5 * relative_limit)).passed is True
    with pytest.raises(AnalysisError, match="does not match cardiac output"):
        validate_physiology_table(_physiology(muscle_flow=1.0 + 1.5 * relative_limit))

    small = 1e-6
    absolute_limit = max(relative * small, absolute)
    assert validate_physiology_table(_flow_table(small, 0.5 * small, 0.5 * small + 0.5 * absolute_limit)).passed
    with pytest.raises(AnalysisError, match="does not match cardiac output"):
        validate_physiology_table(_flow_table(small, 0.5 * small, 0.5 * small + 1.5 * absolute_limit))
    with pytest.raises(ValidationError):
        TissuePhysiology.model_validate(
            {
                "tissue_name": "liver",
                "volume": 1,
                "volume_unit": "L",
                "blood_flow": 1,
                "blood_flow_unit": "L/h",
                "partition_unit": "1",
                "source": "caller note",
            }
        )


def test_missing_liver_negative_baseline_and_relabelled_population_rows() -> None:
    without_liver = _physiology()
    without_liver.tissues[0] = without_liver.tissues[0].model_copy(update={"tissue_name": "kidney"})
    with pytest.raises(AnalysisError, match="liver"):
        steady_state_one_tissue(
            without_liver,
            input_rate=1,
            input_rate_unit="mg/h",
            hepatic_intrinsic_clearance=1,
            clearance_unit="L/h",
            clearance_name="CL",
        )
    with pytest.raises(AnalysisError, match="R0"):
        simulate_linear_turnover(kin=2, kout=1, r0=-1, times_h=[0, 1], context_of_use="research")
    with pytest.raises(AnalysisError, match="relabels evidence grade"):
        load_population_result_table([_population_row(source_evidence_grade="predicted")])
    with pytest.raises(AnalysisError, match="dosing recommendation"):
        load_population_result_table([_population_row(constitutes_dosing_recommendation=False)])


def test_in_process_engines_and_copyleft_imports_stay_refused() -> None:
    with pytest.raises(PlannedCapabilityError):
        run_population()
    with pytest.raises(PlannedCapabilityError):
        run_pbpk()
    with pytest.raises(PlannedCapabilityError):
        run_qsp()
    forbidden = {"nlmixr2", "nlmixr", "rpy2", "pksim", "nonmem", "mobi"}
    root = Path(__file__).resolve().parents[2] / "src" / "pharmacoweave"
    for name in ("population_pk", "pbpk", "qsp"):
        for path in (root / name).rglob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            imported: set[str] = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0].lower() for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0].lower())
            assert imported.isdisjoint(forbidden)


def _estimate(
    subject_id: str,
    name: str,
    value: float,
    *,
    unit: str | None = None,
    standard_error: float | None = None,
) -> IndividualParameterEstimate:
    if unit is None:
        unit = "L/h" if name in {"CL", "CL/F"} else "L"
    return IndividualParameterEstimate(
        subject_id=subject_id,
        parameter_name=name,
        estimate=value,
        unit=unit,
        standard_error=standard_error,
    )


def _population_row(**extra: object) -> dict[str, object]:
    row: dict[str, object] = {
        "subject_id": "SUBJ-01",
        "parameter_name": "CL",
        "estimate": 0.2,
        "unit": "L/h",
        "estimation_method": "caller",
        "evidence_grade": "estimated",
    }
    row.update(extra)
    return row


def _flow_table(cardiac: float, liver_flow: float, muscle_flow: float) -> PhysiologyTable:
    return PhysiologyTable(
        body_weight_kg=0.25,
        body_weight_source="caller body-weight note",
        cardiac_output=cardiac,
        cardiac_output_unit="L/h",
        cardiac_output_source="caller cardiac-output note",
        tissues=[_tissue("liver", liver_flow, 1), _tissue("muscle", muscle_flow, 2)],
    )


def _tissue(name: str, flow: float, partition: float) -> TissuePhysiology:
    return TissuePhysiology(
        tissue_name=name,
        volume=1,
        volume_unit="L",
        blood_flow=flow,
        blood_flow_unit="L/h",
        partition_coefficient=partition,
        partition_unit="1",
        source=f"caller {name} note",
    )


def _physiology(muscle_flow: float = 1) -> PhysiologyTable:
    return PhysiologyTable(
        body_weight_kg=0.25,
        body_weight_source="caller body-weight note",
        cardiac_output=_PBPK["cardiac_output_l_h"],
        cardiac_output_unit="L/h",
        cardiac_output_source="caller cardiac-output note",
        tissues=[
            _tissue("liver", _PBPK["liver_blood_flow_l_h"], _PBPK["kp_liver"]),
            _tissue("muscle", muscle_flow, _PBPK["kp_tissue"]),
        ],
    )
