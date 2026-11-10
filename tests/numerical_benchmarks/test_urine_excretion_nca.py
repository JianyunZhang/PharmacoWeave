"""Urine excretion uses recorded void volume. It does not reuse plasma AUC."""

import json
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from pharmacoweave.cli import main
from pharmacoweave.core.entities import (
    Analyte,
    AnalyteRole,
    Assay,
    Biosample,
    Compound,
    DataClass,
    DoseEvent,
    Observation,
    ObservationStatus,
    ObservationType,
    Route,
    StudyMetadata,
    StudyType,
    Subject,
    Timebase,
    UrineCollection,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import AnalysisError, IngestError, ValidationFailed
from pharmacoweave.ingestion.tabular import write_study_directory
from pharmacoweave.integrations.guards import _PK_MEASUREMENTS
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_nca.runner import load_nca_config, run_nca
from pharmacoweave.pk_nca.urine import PlasmaAucInput, load_urine_nca_config, run_urine_excretion

_REPO = Path(__file__).resolve().parents[2]
_REFERENCE = json.loads((_REPO / "tests" / "regression" / "urine_excretion_quiet.json").read_text(encoding="utf-8"))
_START = datetime(2024, 1, 15, 9, 0, tzinfo=UTC)


def test_source_does_not_reuse_the_plasma_trapezoid() -> None:
    source = (_REPO / "src" / "pharmacoweave" / "pk_nca" / "urine.py").read_text(encoding="utf-8")
    assert "linear_up_log_down" not in source
    assert "auc_last" not in source
    assert "pk_nca.calculate" not in source
    assert "amount_excreted" in _PK_MEASUREMENTS
    assert "renal_clearance" in _PK_MEASUREMENTS


def test_quiet_reference_matches_hand_arithmetic(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    study = _urine_study(
        [
            {
                "start_h": item["interval_start_h"],
                "end_h": item["interval_end_h"],
                "concentration": item["urine_concentration_mg_l"],
                "volume_l": item["urine_volume_l"],
            }
            for item in _REFERENCE["intervals"]
        ]
    )
    auc = PlasmaAucInput(
        source_id="AUC-1",
        subject_id="SUBJ-01",
        analyte_id="ANALYTE-1",
        auc=_REFERENCE["plasma_auc_mg_h_l"],
        auc_unit="mg*h/L",
        matrix="plasma",
    )
    results = run_urine_excretion(study, _REPO / "configs" / "urine_excretion.yaml", plasma_aucs=[auc])
    assert len(results) == 1
    row = results[0]
    assert row.matrix == "urine"
    assert row.amount_excreted_mg == pytest.approx(_REFERENCE["amount_excreted_mg"], rel=1e-12, abs=1e-12)
    assert row.renal_clearance_l_h == pytest.approx(_REFERENCE["renal_clearance_l_h"], rel=1e-12, abs=1e-12)
    assert row.plasma_auc_source_id == "AUC-1"
    assert not hasattr(row, "auc_inf")
    assert not hasattr(row, "cmax")
    assert not hasattr(row, "lambda_z")
    assert any("not a plasma concentration" in item for item in row.limitations)
    for expected, actual in zip(_REFERENCE["intervals"], row.intervals, strict=True):
        assert actual.amount_excreted_mg == pytest.approx(expected["amount_excreted_mg"], rel=1e-12, abs=1e-12)

    assert run_nca(study, _plasma_config()) == []
    bundle = tmp_path / "urine-bundle"
    study.export_bundle(bundle, run_id="urine-quiet", urine_results=results)
    table = pq.read_table(bundle / "analysis" / "urine_excretion.parquet")
    assert "auc_inf" not in table.column_names
    assert "urine_concentration_mg_l" in table.column_names
    assert not (bundle / "analysis" / "nca_parameters.parquet").exists()
    evidence = json.loads((bundle / "evidence" / "pharmacology_evidence.json").read_text(encoding="utf-8"))
    measurements = {item["measurement_type"] for item in evidence["records"]}
    assert measurements <= {"amount_excreted", "renal_clearance"}
    assert "clearance" not in measurements

    directory = tmp_path / "study"
    write_study_directory(study, directory)
    auc_path = tmp_path / "plasma_auc.yaml"
    auc_path.write_text(
        "source_id: AUC-1\nsubject_id: SUBJ-01\nanalyte_id: ANALYTE-1\nauc: 9\nauc_unit: mg*h/L\nmatrix: plasma\n",
        encoding="utf-8",
    )
    assert (
        main(
            [
                "urine-nca",
                str(directory),
                "--config",
                str(_REPO / "configs" / "urine_excretion.yaml"),
                "--plasma-auc",
                str(auc_path),
            ]
        )
        == 0
    )
    printed = capsys.readouterr().out
    assert "amount_excreted_mg=4.5" in printed
    assert "renal_clearance_l_h=0.5" in printed
    with pytest.raises(IngestError):
        load_nca_config(_REPO / "configs" / "urine_excretion.yaml")


def test_incomplete_intervals_keep_known_amounts_and_withhold_clearance() -> None:
    blq = _urine_study(
        [
            {"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25},
            {"start_h": 4, "end_h": 8, "status": ObservationStatus.BLQ, "volume_l": 0.5},
        ]
    )
    blq_row = run_urine_excretion(blq, _config(), plasma_aucs=[_auc()])[0]
    assert blq_row.intervals[0].amount_excreted_mg == pytest.approx(2.5)
    assert blq_row.intervals[1].amount_excreted_mg is None
    assert "not replaced with zero" in (blq_row.intervals[1].omitted_reason or "")
    assert blq_row.amount_excreted_mg is None
    assert blq_row.renal_clearance_l_h is None

    missing = _urine_study(
        [
            {"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25},
            {"start_h": 4, "end_h": 8, "concentration": 4, "volume_l": None},
        ]
    )
    missing_row = run_urine_excretion(missing, _config(), plasma_aucs=[_auc()])[0]
    assert missing_row.intervals[1].urine_concentration_mg_l == pytest.approx(4)
    assert missing_row.intervals[1].amount_excreted_mg is None
    assert "not replaced with zero" in (missing_row.intervals[1].omitted_reason or "")
    assert missing_row.amount_excreted_mg is None
    assert missing_row.renal_clearance_l_h is None

    gap = _urine_study(
        [
            {"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25},
            {"start_h": 6, "end_h": 8, "concentration": 4, "volume_l": 0.5},
        ]
    )
    gap_row = run_urine_excretion(gap, _config(), plasma_aucs=[_auc()])[0]
    assert [item.amount_excreted_mg for item in gap_row.intervals] == pytest.approx([2.5, 2.0])
    assert gap_row.amount_excreted_mg is None
    assert gap_row.renal_clearance_l_h is None

    censored = _urine_study(
        [
            {"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25},
            {"start_h": 4, "end_h": 8, "status": ObservationStatus.ABOVE_ULOQ, "volume_l": 0.5},
        ]
    )
    censored_row = run_urine_excretion(censored, _config(), plasma_aucs=[_auc()])[0]
    assert censored_row.intervals[1].urine_concentration_mg_l is None
    assert censored_row.intervals[1].amount_excreted_mg is None
    assert censored_row.renal_clearance_l_h is None


def test_overlap_wrong_matrix_and_predose_fail() -> None:
    overlap = _urine_study(
        [
            {"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25},
            {"start_h": 2, "end_h": 6, "concentration": 4, "volume_l": 0.5},
        ]
    )
    with pytest.raises(AnalysisError, match="overlapping"):
        run_urine_excretion(overlap, _config(), plasma_aucs=[_auc()])

    plasma = _urine_study([{"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25}], matrix="plasma")
    with pytest.raises(ValidationFailed):
        run_urine_excretion(plasma, _config())
    assert any(item.code == "urine_collection_matrix" for item in validate_study(plasma).issues)

    predose = _urine_study([{"start_h": 0, "end_h": 4, "concentration": 10, "volume_l": 0.25}], predose=True)
    assert any(item.code == "urine_collection_predose" for item in validate_study(predose).issues)
    with pytest.raises(ValueError, match="positive"):
        UrineCollection(
            collection_id="UC-0",
            sample_id="S-0",
            interval_start_since_dose=0,
            interval_start_unit="h",
            volume=0,
            volume_unit="L",
        )


def _config():
    return load_urine_nca_config(_REPO / "configs" / "urine_excretion.yaml")


def _auc() -> PlasmaAucInput:
    return PlasmaAucInput(
        source_id="AUC-1",
        subject_id="SUBJ-01",
        analyte_id="ANALYTE-1",
        auc=9,
        auc_unit="mg*h/L",
        matrix="plasma",
    )


def _plasma_config() -> NcaConfig:
    return NcaConfig(
        analysis="nca",
        design_route="iv_bolus",
        dosing="single_dose",
        sampling="rich_serial",
        auc_method="linear_up_log_down",
        lambda_z_min_points=3,
        lambda_z_min_adjusted_r2=0.9,
        auc_inf_max_extrapolated_fraction=0.2,
        blq_rule="keep_as_censored",
    )


def _urine_study(
    intervals: list[dict[str, object]],
    *,
    matrix: str = "urine",
    predose: bool = False,
) -> PharmacologyStudy:
    metadata = StudyMetadata(
        study_id="PW-URINE",
        protocol_id="SYNTH-URINE",
        title="Synthetic urine excretion",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Rattus norvegicus",
        design_version="1",
        description="Synthetic urine collections. Not an animal experiment.",
        timebases=[Timebase(timebase_id="STUDY-CLOCK", description="Synthetic UTC clock.")],
    )
    samples: list[Biosample] = []
    observations: list[Observation] = []
    collections: list[UrineCollection] = []
    rows = intervals
    if predose:
        rows = [{"start_h": 0, "end_h": -1, "concentration": 10, "volume_l": 0.25, "predose": True}]
    for index, item in enumerate(rows, start=1):
        end_h = float(item["end_h"])  # type: ignore[arg-type]
        sample_id = f"SAMPLE-{index}"
        is_predose = bool(item.get("predose", False))
        samples.append(
            Biosample(
                sample_id=sample_id,
                subject_id="SUBJ-01",
                collection_time_absolute=_START + timedelta(hours=end_h),
                time_since_dose=end_h,
                time_since_dose_unit="h",
                dose_event_id="DOSE-1",
                matrix=matrix,
                timebase_id="STUDY-CLOCK",
                is_predose=is_predose,
            )
        )
        status = item.get("status", ObservationStatus.QUANTIFIED)
        observations.append(
            Observation(
                observation_id=f"OBS-{index}",
                sample_id=sample_id,
                observation_type=ObservationType.CONCENTRATION,
                status=status,  # type: ignore[arg-type]
                analyte_id="ANALYTE-1",
                assay_id="ASSAY-1",
                numeric_value=(
                    None if status is not ObservationStatus.QUANTIFIED else float(item["concentration"])  # type: ignore[arg-type]
                ),
                unit=None if status is not ObservationStatus.QUANTIFIED else "mg/L",
            )
        )
        volume = item.get("volume_l", 0.25)
        collections.append(
            UrineCollection(
                collection_id=f"UC-{index}",
                sample_id=sample_id,
                interval_start_since_dose=float(item["start_h"]),  # type: ignore[arg-type]
                interval_start_unit="h",
                volume=None if volume is None else float(volume),  # type: ignore[arg-type]
                volume_unit=None if volume is None else "L",
            )
        )
    return PharmacologyStudy(
        metadata=metadata,
        subjects=[Subject(subject_id="SUBJ-01", study_id="PW-URINE", species="Rattus norvegicus")],
        compounds=[Compound(compound_id="CMPD-1", name="synthetic")],
        analytes=[
            Analyte(
                analyte_id="ANALYTE-1",
                name="parent",
                analyte_role=AnalyteRole.PARENT,
                matrix=matrix,
                compound_id="CMPD-1",
            )
        ],
        assays=[
            Assay(
                assay_id="ASSAY-1",
                method="synthetic",
                analyte_id="ANALYTE-1",
                matrix=matrix,
                lloq=0.01,
                uloq=1000,
                limit_unit="mg/L",
            )
        ],
        dose_events=[
            DoseEvent(
                dose_event_id="DOSE-1",
                subject_id="SUBJ-01",
                compound_id="CMPD-1",
                absolute_time=_START,
                dose_value=10,
                dose_unit="mg",
                route=Route.IV_BOLUS,
            )
        ],
        biosamples=samples,
        urine_collections=collections,
        observations=observations,
    )
