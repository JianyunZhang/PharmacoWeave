"""Relational validation beyond single-record rules."""

from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import DataClass, ObservationStatus
from pharmacoweave.core.evidence import IdentityLink, ProjectName, ReviewStatus


def test_clean_synthetic_study_passes_strict_validation() -> None:
    report = build_study().check(mode="strict")
    assert report.ok
    assert report.issues == []


def test_missing_subject_is_a_foreign_key_error() -> None:
    study = build_study().model_copy(update={"subjects": []})
    report = study.check()
    assert not report.ok
    assert any(issue.code == "missing_foreign_key" for issue in report.issues)


def test_normalized_dose_must_match_body_weight() -> None:
    study = build_study()
    dose = study.dose_events[0].model_copy(update={"body_weight_normalized_value": 2.0})
    edited = study.model_copy(update={"dose_events": [dose, *study.dose_events[1:]]})
    report = edited.check()
    assert any(issue.code == "dose_body_weight_mismatch" for issue in report.issues)


def test_sample_clock_must_match_the_dose_clock() -> None:
    study = build_study()
    sample = study.biosamples[1].model_copy(update={"time_since_dose": 50.0})
    edited = study.model_copy(update={"biosamples": [study.biosamples[0], sample, *study.biosamples[2:]]})
    report = edited.check()
    assert any(issue.code == "sample_time_inconsistent" for issue in report.issues)


def test_quantified_value_below_lloq_is_rejected() -> None:
    study = build_study()
    observation = study.observations[1].model_copy(update={"numeric_value": 0.001, "raw_value": "0.001"})
    edited = study.model_copy(update={"observations": [study.observations[0], observation, *study.observations[2:]]})
    report = edited.check()
    assert any(issue.code == "below_lloq" for issue in report.issues)


def test_matrix_conflict_and_unlabeled_replicates() -> None:
    study = build_study()
    sample = study.biosamples[1].model_copy(update={"matrix": "brain"})
    duplicate = study.observations[1].model_copy(update={"observation_id": "OBS-DUP"})
    edited = study.model_copy(
        update={
            "biosamples": [study.biosamples[0], sample, *study.biosamples[2:]],
            "observations": [*study.observations, duplicate],
        }
    )
    report = edited.check()
    codes = {issue.code for issue in report.issues}
    assert "matrix_conflict" in codes
    assert "replicate_unlabeled" in codes


def test_unreviewed_identity_is_not_a_match() -> None:
    study = build_study()
    link = IdentityLink(
        link_id="LINK-UNREVIEWED",
        source_project=ProjectName.THERAWEAVE,
        source_study_id="TW-1",
        source_subject_id="SUBJ-01",
        canonical_subject_id="SUBJ-01",
        match_method="text_similarity",
        review_status=ReviewStatus.UNREVIEWED,
        notes="The source id text matches a local id and is still unreviewed.",
    )
    edited = study.model_copy(update={"identity_links": [*study.identity_links, link]})
    research = edited.check(mode="research")
    assert research.ok
    assert any(issue.code == "identity_unreviewed" for issue in research.issues)
    strict = edited.check(mode="strict")
    assert not strict.ok


def test_clinical_records_remain_research_use() -> None:
    study = build_study()
    metadata = study.metadata.model_copy(update={"data_class": DataClass.CLINICAL})
    report = study.model_copy(update={"metadata": metadata}).check()
    assert report.ok
    assert any(issue.code == "clinical_research_use" for issue in report.issues)


def test_blq_rows_in_the_fixture_have_no_number() -> None:
    blq = [item for item in build_study().observations if item.status is ObservationStatus.BLQ]
    assert len(blq) == 4
    assert all(item.numeric_value is None for item in blq)
