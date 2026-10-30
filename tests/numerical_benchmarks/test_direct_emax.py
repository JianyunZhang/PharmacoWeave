"""Noise-free direct Emax recovery. Hill stays fixed at 1."""

import pytest
from examples.direct_effect_and_evidence import E0, EC50, EMAX, build_effects, build_exposures

from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.pkpd import fit_linked_direct_emax

_RELATIVE = 1e-8


def test_direct_emax_recovers_emax_and_ec50() -> None:
    _linked, fit = fit_linked_direct_emax(build_effects(), build_exposures(), e0=E0, e0_unit="percent")
    assert fit.diagnostic.status is DiagnosticStatus.PASS
    found = {item.name: item for item in fit.estimates}
    assert found["Emax"].value == pytest.approx(EMAX, rel=_RELATIVE)
    assert found["Emax"].fixed is False
    assert found["Emax"].unit == "percent"
    assert found["EC50"].value == pytest.approx(EC50, rel=_RELATIVE)
    assert found["EC50"].unit == "mg/L"
    assert found["E0"].fixed is True
    assert found["hill"].value == 1.0
    assert found["hill"].fixed is True
    assert all(point.evidence_grade == "predicted" for point in fit.predicted)
    assert len(fit.predicted) == 5
