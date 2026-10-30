"""The committed synthetic study matches the analytical case."""

import math

from examples.iv_bolus_case import (
    EXAMPLE_DIR,
    LLOQ_UG_PER_ML,
    POSTDOSE_MINUTES,
    SUBJECTS,
    build_study,
    expected_concentration_ug_ml,
)

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.ingestion.tabular import observations_frame, read_study_directory, write_study_directory


def test_committed_example_matches_the_case_and_validates() -> None:
    loaded = PharmacologyStudy.from_files(EXAMPLE_DIR, mode="strict")
    built = build_study()
    assert [item.observation_id for item in loaded.observations] == [item.observation_id for item in built.observations]
    report = loaded.check(mode="strict")
    assert report.ok and report.issues == []
    assert loaded.source_files


def test_postdose_points_follow_the_iv_bolus_equation() -> None:
    study = PharmacologyStudy.from_files(EXAMPLE_DIR)
    by_id = {item.observation_id: item for item in study.observations}
    for subject_id, _sex, body_kg in SUBJECTS:
        concentrations = []
        for minute in POSTDOSE_MINUTES:
            observation = by_id[f"OBS-{subject_id}-T{minute:04d}"]
            expected = expected_concentration_ug_ml(body_kg, minute)
            if expected < LLOQ_UG_PER_ML:
                assert observation.status is ObservationStatus.BLQ
                assert observation.numeric_value is None
            else:
                assert observation.status is ObservationStatus.QUANTIFIED
                assert observation.numeric_value is not None
                assert math.isclose(observation.numeric_value, expected, abs_tol=5e-7)
                concentrations.append(observation.numeric_value)
        assert concentrations


def test_subjects_with_the_same_per_kilogram_parameters_share_concentrations() -> None:
    study = PharmacologyStudy.from_files(EXAMPLE_DIR)
    by_id = {item.observation_id: item for item in study.observations}
    for minute in POSTDOSE_MINUTES:
        left = by_id[f"OBS-SUBJ-01-T{minute:04d}"]
        right = by_id[f"OBS-SUBJ-02-T{minute:04d}"]
        assert left.status is right.status
        assert left.numeric_value == right.numeric_value


def test_directory_roundtrip_preserves_the_study(tmp_path) -> None:
    original = build_study()
    for file_format in ("csv", "tsv", "parquet"):
        directory = tmp_path / file_format
        write_study_directory(original, directory, file_format=file_format)
        loaded = read_study_directory(directory)
        assert [item.subject_id for item in loaded.subjects] == [item.subject_id for item in original.subjects]
        assert [item.numeric_value for item in loaded.observations] == [
            item.numeric_value for item in original.observations
        ]
        assert loaded.check(mode="strict").ok


def test_observations_frame_is_flat() -> None:
    frame = observations_frame(build_study())
    assert len(frame) == 20
    assert "numeric_value" in frame.columns
    assert frame.loc[frame["status"] == "BLQ", "numeric_value"].isna().all()
