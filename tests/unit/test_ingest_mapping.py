"""Unit tests for declarative ingest mapping from raw laboratory files."""

import hashlib
from pathlib import Path

import pytest

from pharmacoweave.core.entities import ObservationStatus, Route
from pharmacoweave.errors import IngestError
from pharmacoweave.ingestion.mapping import (
    IngestMappingConfig,
    ingest_study_from_mapping,
    load_ingest_mapping,
)


def _setup_raw_data(tmp_path: Path) -> Path:
    raw_dir = tmp_path / "raw_lab"
    raw_dir.mkdir(parents=True)

    # 1. Raw subjects file
    (raw_dir / "lab_animals.csv").write_text(
        "AnimalTag,BodyWeightKg,SexCode\n"
        "RAT-01,0.25,F\n"
        "RAT-02,0.30,M\n",
        encoding="utf-8",
    )

    # 2. Raw dosing file
    (raw_dir / "lab_dosing.csv").write_text(
        "DoseTag,AnimalTag,DrugCode,RouteStr,AmountMg,TimeAdmin\n"
        "D-01,RAT-01,TEST-DRUG,iv_bolus,0.25,2024-01-15T09:00:00Z\n"
        "D-02,RAT-02,TEST-DRUG,iv_bolus,0.30,2024-01-15T09:00:00Z\n",
        encoding="utf-8",
    )

    # 3. Raw PK samples and concentrations file
    (raw_dir / "lab_pk_results.csv").write_text(
        "SampleTag,AnimalTag,DoseTag,TimeMin,ConcentrationUgMl,StatusFlag,CollectionTime\n"
        "S-01-01,RAT-01,D-01,5.0,0.95,QUANTIFIED,2024-01-15T09:05:00Z\n"
        "S-01-02,RAT-01,D-01,60.0,0.12,QUANTIFIED,2024-01-15T10:00:00Z\n"
        "S-01-03,RAT-01,D-01,1440.0,,<LLOQ,2024-01-16T09:00:00Z\n"
        "S-02-01,RAT-02,D-02,5.0,1.10,QUANTIFIED,2024-01-15T09:05:00Z\n",
        encoding="utf-8",
    )

    # 4. Ingest mapping YAML
    mapping_yaml = """
schema_name: "pharmacology-ingest-mapping"
schema_version: "0.1.0"
study_metadata:
  study_id: "STUDY-MAPPING-TEST"
  title: "Heterogeneous Raw Data Ingestion Test"
  study_type: "pk"
  data_class: "preclinical_animal"
  species: "Rattus norvegicus"
  design_version: "1.0"
  timebases:
    - timebase_id: "LAB-CLOCK"
      description: "Master NTP clock"

tables:
  subjects:
    source_file: "lab_animals.csv"
    column_mapping:
      AnimalTag: "subject_id"
      BodyWeightKg: "body_weight"
      SexCode: "sex"
    constants:
      species: "Rattus norvegicus"
      body_weight_unit: "kg"
      study_id: "STUDY-MAPPING-TEST"
    value_transforms:
      sex:
        "F": "female"
        "M": "male"

  compounds:
    source_file: "lab_dosing.csv"
    column_mapping:
      DrugCode: "compound_id"
    constants:
      name: "Test Compound"

  analytes:
    source_file: "lab_dosing.csv"
    column_mapping:
      DrugCode: "compound_id"
    constants:
      analyte_id: "ANALYTE-TEST"
      name: "Parent Test Analyte"
      analyte_role: "parent"
      matrix: "plasma"

  assays:
    source_file: "lab_dosing.csv"
    column_mapping: {}
    constants:
      assay_id: "ASSAY-TEST"
      method: "LC-MS/MS"
      analyte_id: "ANALYTE-TEST"
      matrix: "plasma"
      calibration_id: "CAL-01"
      lloq: 0.01
      uloq: 10.0
      limit_unit: "ug/mL"
      qc_status: "accepted"

  dosing_events:
    source_file: "lab_dosing.csv"
    column_mapping:
      DoseTag: "dose_event_id"
      AnimalTag: "subject_id"
      DrugCode: "compound_id"
      RouteStr: "route"
      AmountMg: "dose_value"
      TimeAdmin: "absolute_time"
    constants:
      dose_unit: "mg"

  biosamples:
    source_file: "lab_pk_results.csv"
    column_mapping:
      SampleTag: "sample_id"
      AnimalTag: "subject_id"
      DoseTag: "dose_event_id"
      TimeMin: "time_since_dose"
      CollectionTime: "collection_time_absolute"
    constants:
      time_since_dose_unit: "min"
      matrix: "plasma"
      timebase_id: "LAB-CLOCK"
      is_predose: false

  observations:
    source_file: "lab_pk_results.csv"
    column_mapping:
      SampleTag: "sample_id"
      ConcentrationUgMl: "numeric_value"
      StatusFlag: "status"
    constants:
      observation_type: "concentration"
      unit: "ug/mL"
      analyte_id: "ANALYTE-TEST"
      assay_id: "ASSAY-TEST"
    value_transforms:
      status:
        "QUANTIFIED": "quantified"
        "<LLOQ": "BLQ"
"""
    mapping_path = raw_dir / "ingest_mapping.yaml"
    mapping_path.write_text(mapping_yaml, encoding="utf-8")
    return mapping_path


def test_load_ingest_mapping(tmp_path: Path) -> None:
    mapping_path = _setup_raw_data(tmp_path)
    config = load_ingest_mapping(mapping_path)
    assert isinstance(config, IngestMappingConfig)
    assert config.study_metadata.study_id == "STUDY-MAPPING-TEST"
    assert "subjects" in config.tables
    assert config.tables["subjects"].source_file == "lab_animals.csv"


def test_ingest_study_from_mapping(tmp_path: Path) -> None:
    mapping_path = _setup_raw_data(tmp_path)
    study = ingest_study_from_mapping(mapping_path)

    # 1. Check entities
    assert len(study.subjects) == 2
    assert study.subjects[0].subject_id == "RAT-01"
    assert study.subjects[0].sex == "female"
    assert study.subjects[1].sex == "male"
    assert study.subjects[0].body_weight.value == 0.25
    assert study.subjects[0].body_weight.unit == "kg"

    assert len(study.dose_events) == 2
    assert study.dose_events[0].route == Route.IV_BOLUS
    assert study.dose_events[0].dose_value == 0.25

    assert len(study.biosamples) == 4
    assert len(study.observations) == 4

    # Check value transform for BLQ
    blq_obs = next(o for o in study.observations if o.status == ObservationStatus.BLQ)
    assert blq_obs.numeric_value is None

    # 2. Check source file hashes preserved
    assert len(study.source_files) == 3  # lab_animals.csv, lab_dosing.csv, lab_pk_results.csv
    for sf in study.source_files:
        disk_file = mapping_path.parent / sf.relative_path
        expected_sha = hashlib.sha256(disk_file.read_bytes()).hexdigest()
        assert sf.sha256 == expected_sha

    # 3. Check study check / validation
    report = study.check()
    assert report.ok, str(report)


def test_ingest_study_missing_file_raises(tmp_path: Path) -> None:
    mapping_path = _setup_raw_data(tmp_path)
    # Delete one file
    (mapping_path.parent / "lab_animals.csv").unlink()

    with pytest.raises(IngestError, match="not found"):
        ingest_study_from_mapping(mapping_path)
