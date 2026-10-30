"""Toxicology records stay observations."""

import pytest
from pydantic import ValidationError

from pharmacoweave.core.evidence import TOX_RESEARCH_LIMITATION
from pharmacoweave.errors import ContractViolation
from pharmacoweave.toxicology import ToxEndpoint, ToxKind, ToxSeverity, dosing_recommendation


def _endpoint(**overrides: object) -> ToxEndpoint:
    payload: dict[str, object] = {
        "endpoint_id": "TOX-1",
        "compound_id": "PW-0001",
        "kind": ToxKind.HISTOPATHOLOGY,
        "finding": "Minimal hypertrophy in a synthetic example.",
        "severity": ToxSeverity.MINIMAL,
        "limitations": [TOX_RESEARCH_LIMITATION],
        "dose_event_id": "DOSE-1",
    }
    payload.update(overrides)
    return ToxEndpoint(**payload)  # type: ignore[arg-type]


def test_endpoint_links_to_dose_or_exposure_and_grades_histopathology() -> None:
    assert _endpoint().evidence_grade == "observed"
    linked = _endpoint(dose_event_id=None, exposure_ref="nca:SUBJ-01:auc_inf")
    assert linked.exposure_ref == "nca:SUBJ-01:auc_inf"
    with pytest.raises(ValidationError, match="dose_event_id or an exposure_ref"):
        _endpoint(dose_event_id=None, exposure_ref=None)
    with pytest.raises(ValidationError, match="severity"):
        _endpoint(severity=ToxSeverity.NOT_GRADED)


def test_incidence_dose_units_and_clinical_advice_are_rejected() -> None:
    with pytest.raises(ValidationError, match="incidence_affected cannot exceed"):
        _endpoint(
            kind=ToxKind.MORTALITY,
            severity=ToxSeverity.MARKED,
            incidence_affected=3,
            incidence_examined=2,
        )
    mortality = _endpoint(
        kind=ToxKind.MORTALITY,
        severity=ToxSeverity.MARKED,
        incidence_affected=1,
        incidence_examined=4,
    )
    assert mortality.incidence_affected == 1
    with pytest.raises(ValidationError, match="administered dose"):
        _endpoint(
            kind=ToxKind.BODY_WEIGHT_CHANGE,
            severity=ToxSeverity.NOT_GRADED,
            numeric_value=-2.0,
            unit="mg",
        )
    change = _endpoint(
        kind=ToxKind.BODY_WEIGHT_CHANGE,
        severity=ToxSeverity.NOT_GRADED,
        numeric_value=-4.0,
        unit="percent",
    )
    assert change.numeric_value == -4.0
    with pytest.raises(ValidationError):
        _endpoint(recommendation="10 mg")  # type: ignore[call-arg]
    with pytest.raises(ContractViolation, match="clinical dose"):
        dosing_recommendation(mortality)
