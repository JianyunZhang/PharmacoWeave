"""Original-scale damped nonlinear least squares for direct and sigmoid Emax."""

import pytest
from examples.noisy_emax_case import (
    DIRECT_RECOVERY_RELATIVE_TOLERANCE,
    E0,
    SIGMOID_RECOVERY_RELATIVE_TOLERANCE,
    SIGMOID_TRUE_EC50,
    SIGMOID_TRUE_EMAX,
    SIGMOID_TRUE_HILL,
    TRUE_EC50,
    TRUE_EMAX,
    build_noisy_direct_case,
    build_noisy_sigmoid_case,
)

from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pd_models import (
    DirectEffectPair,
    EffectDirection,
    emax_effect,
    estimate_hill,
    fit_direct_emax,
    fit_sigmoid_emax,
)
from pharmacoweave.pkpd import fit_linked_direct_emax


def _pair(index: int, concentration: float, effect: float) -> DirectEffectPair:
    return DirectEffectPair(
        effect_id=f"EFF-{index}",
        exposure_ref=f"CONC-{index}",
        concentration=concentration,
        concentration_unit="mg/L",
        effect=effect,
        effect_unit="percent",
    )


def test_noisy_direct_emax_recovers_across_baseline_dip_and_reports_ci() -> None:
    effects, exposures = build_noisy_direct_case()
    linked, fit = fit_linked_direct_emax(effects, exposures, e0=E0, e0_unit="percent")
    assert [item.effect_id for item in linked.omitted] == ["NOISY-EFF-BLQ"]
    assert fit.fit.converged is True
    assert fit.diagnostic.status is DiagnosticStatus.WARNING
    found = {item.name: item for item in fit.estimates}
    assert found["Emax"].value == pytest.approx(TRUE_EMAX, rel=DIRECT_RECOVERY_RELATIVE_TOLERANCE)
    assert found["EC50"].value == pytest.approx(TRUE_EC50, rel=DIRECT_RECOVERY_RELATIVE_TOLERANCE)
    assert found["hill"].fixed is True
    assert found["hill"].value == 1.0
    for name, truth in (("Emax", TRUE_EMAX), ("EC50", TRUE_EC50)):
        param = found[name]
        assert param.standard_error is not None and param.standard_error > 0
        assert param.ci_lower is not None and param.ci_upper is not None
        assert param.ci_level == 0.95
        assert param.ci_lower < param.value < param.ci_upper
        assert param.ci_lower <= truth <= param.ci_upper
    assert found["EC50"].ci_lower > 0

    legacy = fit_direct_emax(linked.pairs, e0=E0, e0_unit="percent", method="linearized")
    assert legacy.diagnostic.status is DiagnosticStatus.FAIL
    assert legacy.estimates == ()


@pytest.mark.parametrize("true_hill", [0.75, 2.2])
def test_noise_free_sigmoid_emax_recovers_hill_and_ec50(true_hill: float) -> None:
    concentrations = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
    pairs = [
        _pair(index, conc, emax_effect(concentration=conc, e0=10.0, emax=80.0, ec50=8.0, hill=true_hill))
        for index, conc in enumerate(concentrations, start=1)
    ]
    fit = fit_sigmoid_emax(pairs, e0=10.0, e0_unit="percent", direction=EffectDirection.INCREASE)
    assert fit.fit.converged is True
    assert fit.diagnostic.status is DiagnosticStatus.PASS
    found = {item.name: item for item in fit.estimates}
    assert found["Emax"].value == pytest.approx(80.0, rel=1e-6)
    assert found["EC50"].value == pytest.approx(8.0, rel=1e-6)
    assert found["hill"].value == pytest.approx(true_hill, rel=1e-6)
    assert found["hill"].fixed is False
    assert found["hill"].standard_error == pytest.approx(0.0, abs=1e-9)

    via_alias = estimate_hill(pairs, e0=10.0, e0_unit="percent", direction=EffectDirection.INCREASE)
    alias_found = {item.name: item for item in via_alias.estimates}
    assert alias_found["hill"].value == pytest.approx(true_hill, rel=1e-6)


def test_noisy_sigmoid_emax_recovers_parameters_and_positive_intervals() -> None:
    effects, exposures = build_noisy_sigmoid_case()
    _linked, fit = fit_linked_direct_emax(
        effects,
        exposures,
        e0=E0,
        e0_unit="percent",
        estimate_hill_coefficient=True,
    )
    assert fit.fit.converged is True
    assert fit.diagnostic.status in {DiagnosticStatus.PASS, DiagnosticStatus.WARNING}
    found = {item.name: item for item in fit.estimates}
    assert found["Emax"].value == pytest.approx(SIGMOID_TRUE_EMAX, rel=SIGMOID_RECOVERY_RELATIVE_TOLERANCE)
    assert found["EC50"].value == pytest.approx(SIGMOID_TRUE_EC50, rel=SIGMOID_RECOVERY_RELATIVE_TOLERANCE)
    assert found["hill"].value == pytest.approx(SIGMOID_TRUE_HILL, rel=SIGMOID_RECOVERY_RELATIVE_TOLERANCE)
    for name in ("Emax", "EC50", "hill"):
        param = found[name]
        assert param.standard_error is not None and param.standard_error > 0
        assert param.ci_lower is not None and param.ci_upper is not None
        assert param.ci_lower < param.value < param.ci_upper
    assert found["EC50"].ci_lower > 0
    assert found["hill"].ci_lower > 0


def test_unidentifiable_plateau_and_step_curves_fail_without_ec50_or_hill() -> None:
    linear_only = [
        _pair(index, conc, emax_effect(concentration=conc, e0=10.0, emax=80.0, ec50=100.0, hill=1.0))
        for index, conc in enumerate((0.05, 0.1, 0.2, 0.4), start=1)
    ]
    no_plateau = fit_direct_emax(linear_only, e0=10.0, e0_unit="percent")
    assert no_plateau.diagnostic.status is DiagnosticStatus.FAIL
    assert no_plateau.estimates == ()
    assert no_plateau.predicted == ()

    saturated_only = [
        _pair(index, conc, emax_effect(concentration=conc, e0=10.0, emax=80.0, ec50=0.5, hill=1.0))
        for index, conc in enumerate((200.0, 400.0, 800.0, 1600.0), start=1)
    ]
    already_flat = fit_direct_emax(saturated_only, e0=10.0, e0_unit="percent")
    assert already_flat.diagnostic.status is DiagnosticStatus.FAIL
    assert already_flat.estimates == ()

    step_gap = [
        _pair(1, 0.1, 10.0),
        _pair(2, 0.2, 10.0),
        _pair(3, 100.0, 90.0),
        _pair(4, 200.0, 90.0),
    ]
    degenerate_sigmoid = fit_sigmoid_emax(step_gap, e0=10.0, e0_unit="percent")
    assert degenerate_sigmoid.diagnostic.status is DiagnosticStatus.FAIL
    assert degenerate_sigmoid.estimates == ()

    too_short_sigmoid = [
        _pair(index, conc, emax_effect(concentration=conc, e0=10.0, emax=80.0, ec50=8.0, hill=2.0))
        for index, conc in enumerate((2.0, 8.0, 32.0), start=1)
    ]
    short_sigmoid = fit_sigmoid_emax(too_short_sigmoid, e0=10.0, e0_unit="percent")
    assert short_sigmoid.diagnostic.status is DiagnosticStatus.FAIL
    assert "four" in short_sigmoid.reason.lower()

    with pytest.raises(AnalysisError, match="Linearized"):
        fit_direct_emax(linear_only, e0=10.0, e0_unit="percent", method="linearized", estimate_hill_coefficient=True)
