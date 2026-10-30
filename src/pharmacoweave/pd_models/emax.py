"""Direct-effect and sigmoid Emax estimation with fixed baseline E0.

By default, ``fit_direct_emax`` uses damped Levenberg–Marquardt nonlinear least
squares on the original effect scale:

    E = E0 + Emax * C^hill / (EC50^hill + C^hill)

with ``hill`` fixed at 1 unless ``estimate_hill_coefficient=True`` (or
``fit_sigmoid_emax`` / ``estimate_hill`` is called). The legacy double-reciprocal
rearrangement remains available via ``method='linearized'`` for ``hill=1``.
Non-identifiable curves fail without reporting EC50 or hill. Model predictions
stay labeled ``predicted`` and cannot be stored as observed effects.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, PositiveFloat, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus, Metric, ModelFit, ModelSpec
from pharmacoweave.core.numerics import (
    invert_symmetric,
    normalized_curvature_ratio,
    solve_linear_system,
    student_t_two_sided_quantile,
)
from pharmacoweave.core.units import CONC_AMOUNT, CONC_MASS, convert, parse_unit
from pharmacoweave.errors import AnalysisError, ContractViolation, UnitError
from pharmacoweave.pd_models.records import EffectDirection, EffectObservation
from pharmacoweave.pd_models.specs import DIRECT_EMAX, SIGMOID_EMAX
from pharmacoweave.version import SCHEMA_VERSION, __version__

_METHOD_NLS = "direct_emax_nls_v1"
_METHOD_SIGMOID_NLS = "sigmoid_emax_nls_v1"
_METHOD_LINEARIZED = "direct_emax_linearized_v1"
_CI_LEVEL = 0.95

_NLS_LIMITATIONS_HILL1 = (
    RESEARCH_USE_LIMITATION,
    "Direct Emax with a fixed baseline and hill fixed at 1, estimated on the original effect scale.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton approximation.",
    "Each effect is paired with an explicit concentration. AUC, Cmax, lag, and indirect response are not used.",
    "A failed diagnostic does not report EC50.",
)
_NLS_LIMITATIONS_SIGMOID = (
    RESEARCH_USE_LIMITATION,
    "Sigmoid direct Emax with a fixed baseline and estimated Hill coefficient on the original effect scale.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton approximation.",
    "Each effect is paired with an explicit concentration. AUC, Cmax, lag, and indirect response are not used.",
    "A failed diagnostic does not report EC50 or hill.",
)
_LINEARIZED_LIMITATIONS = (
    RESEARCH_USE_LIMITATION,
    (
        "Direct Emax linearized by double-reciprocal OLS with fixed baseline "
        "and hill fixed at 1. No standard error is reported."
    ),
    "Each effect is paired with an explicit concentration. AUC, Cmax, lag, and indirect response are not used.",
    "A failed diagnostic does not report EC50.",
)


class EffectParameter(Record):
    """One direct-effect parameter. EC50 is a concentration, not a clearance."""

    fit_id: Identifier
    name: str = Field(min_length=1)
    value: FiniteFloat
    unit: str
    fixed: bool
    estimation_method: str = Field(min_length=1)
    standard_error: float | None = Field(default=None, ge=0)
    ci_lower: float | None = None
    ci_upper: float | None = None
    ci_level: float | None = Field(default=None, gt=0, lt=1)

    @model_validator(mode="after")
    def unit_matches_name(self) -> "EffectParameter":
        parsed = parse_unit(self.unit)
        if self.name == "EC50":
            if self.value <= 0:
                raise ValueError("EC50 must be positive.")
            if parsed.dimensions not in {CONC_MASS, CONC_AMOUNT}:
                raise ValueError("EC50 must use a concentration unit.")
            if self.ci_lower is not None and self.ci_lower <= 0:
                raise ValueError("EC50 ci_lower must be positive.")
        if self.name == "hill":
            if self.value <= 0:
                raise ValueError("hill must be positive.")
            if self.unit != "1":
                raise ValueError("hill uses the dimensionless unit 1.")
            if self.ci_lower is not None and self.ci_lower <= 0:
                raise ValueError("hill ci_lower must be positive.")
        if (self.ci_lower is None) != (self.ci_upper is None):
            raise ValueError("ci_lower and ci_upper are recorded together.")
        if self.ci_lower is not None and self.ci_upper is not None and self.ci_lower > self.ci_upper:
            raise ValueError("ci_lower must be less than or equal to ci_upper.")
        if self.ci_lower is not None and self.ci_level is None:
            raise ValueError("A confidence interval requires ci_level.")
        return self


class DirectEffectPair(Record):
    """One observed effect and the concentration it is paired with."""

    effect_id: Identifier
    exposure_ref: Identifier
    concentration: PositiveFloat
    concentration_unit: str
    effect: FiniteFloat
    effect_unit: str

    @model_validator(mode="after")
    def units(self) -> "DirectEffectPair":
        if parse_unit(self.concentration_unit).dimensions not in {CONC_MASS, CONC_AMOUNT}:
            raise ValueError("concentration_unit must be a concentration.")
        parse_unit(self.effect_unit)
        return self


@dataclass(frozen=True)
class PredictedEffect:
    """A model effect. ``evidence_grade`` stays predicted."""

    effect_id: str
    exposure_ref: str
    value: float
    unit: str
    derived_from: tuple[str, ...]
    evidence_grade: str = "predicted"

    def __post_init__(self) -> None:
        if self.evidence_grade != "predicted":
            raise ValueError("A model effect is predicted and cannot be labeled observed.")


class PredictedEffectRecord(Record):
    """Serializable schema contract for a model-predicted effect."""

    fit_id: Identifier
    effect_id: Identifier
    exposure_ref: Identifier
    value: FiniteFloat
    unit: str
    derived_from: list[str] = Field(min_length=1)
    evidence_grade: Literal["predicted"] = "predicted"


class PdAnalysisRecords(Record):
    """Machine-readable PD contract exported to ``schemas/pd_records.schema.json``."""

    effect_observations: list[EffectObservation] = Field(default_factory=list)
    direct_effect_pairs: list[DirectEffectPair] = Field(default_factory=list)
    model_specs: list[ModelSpec] = Field(default_factory=list)
    model_fits: list[ModelFit] = Field(default_factory=list)
    effect_parameters: list[EffectParameter] = Field(default_factory=list)
    predicted_effects: list[PredictedEffectRecord] = Field(default_factory=list)
    diagnostics: list[DiagnosticRecord] = Field(default_factory=list)
    schema_name: Literal["pharmacology-pd-records"] = "pharmacology-pd-records"
    schema_version: str = SCHEMA_VERSION

    @model_validator(mode="after")
    def check_version(self) -> "PdAnalysisRecords":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        return self


@dataclass(frozen=True)
class DirectEmaxFit:
    """One direct or sigmoid Emax fit. Failed fits have empty estimates and predictions."""

    fit: ModelFit
    estimates: tuple[EffectParameter, ...]
    diagnostic: DiagnosticRecord
    predicted: tuple[PredictedEffect, ...]
    effect_ids: tuple[str, ...]
    exposure_refs: tuple[str, ...]
    reason: str
    pairs: tuple[DirectEffectPair, ...] = ()
    e0: float | None = None
    e0_unit: str | None = None


def emax_effect(*, concentration: float, e0: float, emax: float, ec50: float, hill: float = 1.0) -> float:
    """Return E0 + Emax * C^hill / (EC50^hill + C^hill)."""

    if concentration < 0 or ec50 <= 0 or hill <= 0 or not math.isfinite(e0) or not math.isfinite(emax):
        raise ValueError("concentration must be >= 0, EC50 and hill must be > 0, and E0 and Emax must be finite.")
    if concentration == 0:
        return e0
    z = hill * (math.log(concentration) - math.log(ec50))
    return e0 + emax * _logistic(z)


def fit_sigmoid_emax(
    pairs: Sequence[DirectEffectPair],
    *,
    e0: float,
    e0_unit: str,
    direction: EffectDirection = EffectDirection.UNSPECIFIED,
    fit_id: str = "sigmoid-emax",
) -> DirectEmaxFit:
    """Estimate Emax, EC50, and Hill coefficient on the original effect scale."""

    return fit_direct_emax(
        pairs,
        e0=e0,
        e0_unit=e0_unit,
        direction=direction,
        fit_id=fit_id,
        method="nonlinear_ls",
        estimate_hill_coefficient=True,
    )


def estimate_hill(
    pairs: Sequence[DirectEffectPair],
    *,
    e0: float,
    e0_unit: str,
    direction: EffectDirection = EffectDirection.UNSPECIFIED,
    fit_id: str = "sigmoid-emax",
) -> DirectEmaxFit:
    """Estimate Emax, EC50, and the Sigmoid Hill coefficient."""

    return fit_sigmoid_emax(
        pairs,
        e0=e0,
        e0_unit=e0_unit,
        direction=direction,
        fit_id=fit_id,
    )


def effect_observation_from_prediction(predicted: PredictedEffect) -> None:
    """Refuse to store a prediction as an observed effect."""

    raise ContractViolation(
        f"{predicted.effect_id} is a model prediction and cannot be stored as an observed effect."
    )


def fit_direct_emax(
    pairs: Sequence[DirectEffectPair],
    *,
    e0: float,
    e0_unit: str,
    direction: EffectDirection = EffectDirection.UNSPECIFIED,
    fit_id: str = "direct-emax",
    method: Literal["nonlinear_ls", "linearized"] = "nonlinear_ls",
    estimate_hill_coefficient: bool = False,
) -> DirectEmaxFit:
    """Estimate direct or sigmoid Emax. Non-identifiable curves fail without EC50."""

    if method not in {"nonlinear_ls", "linearized"}:
        raise AnalysisError("method must be 'nonlinear_ls' or 'linearized'.")
    if method == "linearized" and estimate_hill_coefficient:
        raise AnalysisError("Linearized Emax fixes hill at 1; use method='nonlinear_ls' to estimate hill.")
    if not math.isfinite(e0):
        raise AnalysisError("e0 must be finite.")
    parse_unit(e0_unit)

    rows = tuple(pairs)
    effect_ids = tuple(row.effect_id for row in rows)
    exposure_refs = tuple(row.exposure_ref for row in rows)
    if len(effect_ids) != len(set(effect_ids)):
        raise AnalysisError("effect_id values must be unique inside one direct Emax fit.")
    if len(exposure_refs) != len(set(exposure_refs)):
        raise AnalysisError("exposure_ref values must be unique inside one direct Emax fit.")

    spec = SIGMOID_EMAX if estimate_hill_coefficient else DIRECT_EMAX
    base_limitations: tuple[str, ...]
    if method == "linearized":
        method_tag = _METHOD_LINEARIZED
        base_limitations = _LINEARIZED_LIMITATIONS
    elif estimate_hill_coefficient:
        method_tag = _METHOD_SIGMOID_NLS
        base_limitations = _NLS_LIMITATIONS_SIGMOID
    else:
        method_tag = _METHOD_NLS
        base_limitations = _NLS_LIMITATIONS_HILL1

    min_points = 4 if estimate_hill_coefficient else 3
    if len(rows) < min_points:
        label = "four" if estimate_hill_coefficient else "three"
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            f"Fewer than {label} quantified concentration-effect pairs.",
        )

    concentration_unit = rows[0].concentration_unit
    try:
        concentrations = [convert(row.concentration, row.concentration_unit, concentration_unit) for row in rows]
        effects = [convert(row.effect, row.effect_unit, e0_unit) for row in rows]
    except UnitError as exc:
        raise AnalysisError(str(exc)) from exc

    normalized_pairs = tuple(
        DirectEffectPair(
            effect_id=row.effect_id,
            exposure_ref=row.exposure_ref,
            concentration=conc,
            concentration_unit=concentration_unit,
            effect=eff,
            effect_unit=e0_unit,
        )
        for row, conc, eff in zip(rows, concentrations, effects, strict=True)
    )

    if method == "linearized":
        return _fit_linearized(
            fit_id=fit_id,
            rows=normalized_pairs,
            concentrations=concentrations,
            effects=effects,
            concentration_unit=concentration_unit,
            e0=e0,
            e0_unit=e0_unit,
            direction=direction,
            effect_ids=effect_ids,
            exposure_refs=exposure_refs,
        )

    return _fit_nonlinear_ls(
        fit_id=fit_id,
        spec=spec,
        method_tag=method_tag,
        base_limitations=base_limitations,
        rows=normalized_pairs,
        concentrations=concentrations,
        effects=effects,
        concentration_unit=concentration_unit,
        e0=e0,
        e0_unit=e0_unit,
        direction=direction,
        estimate_hill_coefficient=estimate_hill_coefficient,
        effect_ids=effect_ids,
        exposure_refs=exposure_refs,
    )


def _fit_nonlinear_ls(
    *,
    fit_id: str,
    spec: ModelSpec,
    method_tag: str,
    base_limitations: tuple[str, ...],
    rows: tuple[DirectEffectPair, ...],
    concentrations: list[float],
    effects: list[float],
    concentration_unit: str,
    e0: float,
    e0_unit: str,
    direction: EffectDirection,
    estimate_hill_coefficient: bool,
    effect_ids: tuple[str, ...],
    exposure_refs: tuple[str, ...],
) -> DirectEmaxFit:
    scale = max(1.0, abs(e0), max(abs(value) for value in effects))
    deltas = [value - e0 for value in effects]
    max_abs_delta = max(abs(delta) for delta in deltas)
    span = max(deltas) - min(deltas)
    if max_abs_delta <= 1e-12 * scale or span <= 1e-12 * max(1.0, max_abs_delta):
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Effects do not vary with concentration, so Emax and EC50 are not identifiable.",
        )

    distinct_c = {round(math.log(c), 9) for c in concentrations}
    min_distinct = 4 if estimate_hill_coefficient else 2
    if len(distinct_c) < min_distinct or max(concentrations) - min(concentrations) <= 0:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Concentrations do not vary across enough distinct levels to identify the Emax curve.",
        )

    if min(deltas) < -0.25 * max_abs_delta and max(deltas) > 0.25 * max_abs_delta:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            (
                "Effects move substantially both above and below the fixed baseline, "
                "so one monotonic Emax is not identifiable."
            ),
        )

    solved = _solve_emax_lm(
        concentrations=concentrations,
        deltas=deltas,
        estimate_hill_coefficient=estimate_hill_coefficient,
    )
    if solved is None:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Nonlinear least squares did not converge or parameters were non-identifiable.",
        )
    emax, ec50, hill, sse, iterations, hessian, curvature = solved

    if ec50 > 1.25 * max(concentrations):
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            (
                "Concentrations do not reach the Emax transition/plateau "
                "(EC50 exceeds tested concentrations), so Emax and EC50 are not separately identifiable."
            ),
        )
    if ec50 < 0.05 * min(concentrations):
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "All concentrations lie far above estimated EC50 on the plateau, so EC50 is not identifiable.",
        )
    if estimate_hill_coefficient:
        if not 0.15 <= hill <= 8.0:
            return _failed(
                fit_id,
                spec,
                method_tag,
                base_limitations,
                effect_ids,
                exposure_refs,
                rows,
                e0,
                e0_unit,
                "Hill coefficient hit the identifiability boundary, indicating a degenerate step or flat curve.",
            )
        in_transition = sum(0.10 * abs(emax) <= abs(delta) <= 0.90 * abs(emax) for delta in deltas)
        if in_transition < 1:
            return _failed(
                fit_id,
                spec,
                method_tag,
                base_limitations,
                effect_ids,
                exposure_refs,
                rows,
                e0,
                e0_unit,
                (
                    "Sigmoid Emax requires at least one observation in the transition zone "
                    "(10%-90% of Emax) to identify hill."
                ),
            )
    if curvature < 1e-5:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "The concentration-effect curvature matrix is ill-conditioned, so EC50 is not identifiable.",
        )

    if direction is EffectDirection.INCREASE and emax <= 0:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "The declared direction is increase, but Emax is not positive. EC50 was not reported.",
        )
    if direction is EffectDirection.DECREASE and emax >= 0:
        return _failed(
            fit_id,
            spec,
            method_tag,
            base_limitations,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "The declared direction is decrease, but Emax is not negative. EC50 was not reported.",
        )

    predicted_values = [
        emax_effect(concentration=value, e0=e0, emax=emax, ec50=ec50, hill=hill) for value in concentrations
    ]
    residual = max(abs(pred - observed) for pred, observed in zip(predicted_values, effects, strict=True))
    status = DiagnosticStatus.PASS if residual <= 1e-4 * scale else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["An effect residual exceeds the noise-free tolerance."]

    n_free = 3 if estimate_hill_coefficient else 2
    df = len(rows) - n_free
    uncertainties = _parameter_uncertainties(
        emax=emax,
        ec50=ec50,
        hill=hill,
        sse=sse,
        df=df,
        hessian=hessian,
        estimate_hill_coefficient=estimate_hill_coefficient,
        warnings=warnings,
    )
    emax_unc, ec50_unc, hill_unc = uncertainties

    estimates = (
        EffectParameter(
            fit_id=fit_id,
            name="E0",
            value=e0,
            unit=e0_unit,
            fixed=True,
            estimation_method="fixed",
        ),
        EffectParameter(
            fit_id=fit_id,
            name="Emax",
            value=emax,
            unit=e0_unit,
            fixed=False,
            estimation_method=method_tag,
            **emax_unc,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="EC50",
            value=ec50,
            unit=concentration_unit,
            fixed=False,
            estimation_method=method_tag,
            **ec50_unc,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="hill",
            value=hill,
            unit="1",
            fixed=not estimate_hill_coefficient,
            estimation_method=method_tag if estimate_hill_coefficient else "fixed",
            **hill_unc,
        ),
    )
    predicted = tuple(
        PredictedEffect(
            effect_id=row.effect_id,
            exposure_ref=row.exposure_ref,
            value=value,
            unit=e0_unit,
            derived_from=(f"fit:{fit_id}", f"exposure:{row.exposure_ref}"),
        )
        for row, value in zip(rows, predicted_values, strict=True)
    )
    metrics = [
        Metric(name="n_points", value=float(len(rows))),
        Metric(name="n_iterations", value=float(iterations)),
        Metric(name="sse", value=sse),
        Metric(name="curvature_ratio", value=curvature),
        Metric(name="max_abs_residual", value=residual, unit=e0_unit),
        Metric(name="Emax", value=emax, unit=e0_unit),
        Metric(name="EC50", value=ec50, unit=concentration_unit),
        Metric(name="hill", value=hill, unit="1"),
    ]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(base_limitations),
        warnings=warnings,
        metrics=metrics,
    )
    fit = _model_fit(fit_id, spec, method_tag, status, converged=True, limitations=list(base_limitations))
    return DirectEmaxFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        predicted=predicted,
        effect_ids=effect_ids,
        exposure_refs=exposure_refs,
        reason="",
        pairs=rows,
        e0=e0,
        e0_unit=e0_unit,
    )


def _solve_emax_lm(
    *,
    concentrations: list[float],
    deltas: list[float],
    estimate_hill_coefficient: bool,
) -> tuple[float, float, float, float, int, list[list[float]], float] | None:
    log_c = [math.log(c) for c in concentrations]
    lo, hi = min(log_c), max(log_c)
    span = max(hi - lo, 0.5)
    eta_grid = [lo - 0.5 * span + (k / 6.0) * (2.0 * span) for k in range(7)]
    xi_grid = [math.log(val) for val in (0.5, 0.8, 1.0, 1.5, 2.2, 3.5)] if estimate_hill_coefficient else [0.0]

    def sse_at(emax_val: float, eta_val: float, xi_val: float) -> float:
        hill_val = math.exp(xi_val)
        total = 0.0
        for lc, delta in zip(log_c, deltas, strict=True):
            pred = emax_val * _logistic(hill_val * (lc - eta_val))
            diff = delta - pred
            total += diff * diff
        return total

    best_start: tuple[float, float, float, float] | None = None
    for eta0 in eta_grid:
        for xi0 in xi_grid:
            hill0 = math.exp(xi0)
            gs = [_logistic(hill0 * (lc - eta0)) for lc in log_c]
            denom = sum(g * g for g in gs)
            if denom <= 1e-18:
                continue
            emax0 = sum(delta * g for delta, g in zip(deltas, gs, strict=True)) / denom
            if abs(emax0) <= 1e-15 or not math.isfinite(emax0):
                continue
            cost = sse_at(emax0, eta0, xi0)
            if math.isfinite(cost) and (best_start is None or cost < best_start[0]):
                best_start = (cost, emax0, eta0, xi0)

    if best_start is None:
        return None
    current, emax, eta, xi = best_start
    p = 3 if estimate_hill_coefficient else 2

    def normal_equations(emax_val: float, eta_val: float, xi_val: float) -> tuple[list[list[float]], list[float]]:
        hill_val = math.exp(xi_val)
        h_mat = [[0.0] * p for _ in range(p)]
        grad = [0.0] * p
        for lc, delta in zip(log_c, deltas, strict=True):
            z = hill_val * (lc - eta_val)
            g = _logistic(z)
            dg = g * (1.0 - g)
            residual = delta - emax_val * g
            col0 = g
            col1 = -emax_val * hill_val * dg
            cols = [col0, col1, emax_val * dg * z] if estimate_hill_coefficient else [col0, col1]
            for i in range(p):
                grad[i] += cols[i] * residual
                for j in range(i, p):
                    h_mat[i][j] += cols[i] * cols[j]
        for i in range(p):
            for j in range(i):
                h_mat[i][j] = h_mat[j][i]
        return h_mat, grad

    lam = 1e-3
    converged = False
    iterations = 0
    for _step in range(80):
        iterations += 1
        hessian, grad = normal_equations(emax, eta, xi)
        scaled_grad = max(
            abs(grad[i]) / max(1.0, math.sqrt(max(0.0, hessian[i][i]))) for i in range(p)
        )
        if scaled_grad <= 1e-10 * (1.0 + math.sqrt(current)):
            converged = True
            break
        accepted = False
        for _attempt in range(20):
            damped = [row[:] for row in hessian]
            for i in range(p):
                scale_i = hessian[i][i] if hessian[i][i] > 1e-12 else 1.0
                damped[i][i] += lam * scale_i
            step = solve_linear_system(damped, grad)
            if step is None:
                lam = min(lam * 8.0, 1e12)
                continue
            step_norm = max(
                abs(step[0]) / max(1.0, abs(emax)),
                *(abs(step[k]) for k in range(1, p)),
            )
            if step_norm <= 1e-10:
                converged = True
                accepted = True
                break
            trial_emax = emax + step[0]
            trial_eta = eta + step[1]
            trial_xi = xi + step[2] if estimate_hill_coefficient else 0.0
            if not math.isfinite(trial_emax) or abs(trial_eta) > 30.0 or abs(trial_xi) > 3.5:
                lam = min(lam * 8.0, 1e12)
                continue
            trial = sse_at(trial_emax, trial_eta, trial_xi)
            if not math.isfinite(trial):
                lam = min(lam * 8.0, 1e12)
                continue
            if trial <= current * (1.0 + 1e-14) + 1e-24:
                emax, eta, xi, current = trial_emax, trial_eta, trial_xi, trial
                lam = max(lam / 3.0, 1e-14)
                accepted = True
                break
            lam = min(lam * 8.0, 1e12)
        if converged or not accepted:
            break

    if not converged or abs(emax) <= 1e-15 or not math.isfinite(emax):
        return None
    ec50 = math.exp(eta)
    hill = math.exp(xi) if estimate_hill_coefficient else 1.0
    if not (math.isfinite(ec50) and ec50 > 0 and math.isfinite(hill) and hill > 0 and math.isfinite(current)):
        return None
    final_hessian, _ = normal_equations(emax, eta, xi)
    curvature = normalized_curvature_ratio(final_hessian)
    return emax, ec50, hill, current, iterations, final_hessian, curvature


def _parameter_uncertainties(
    *,
    emax: float,
    ec50: float,
    hill: float,
    sse: float,
    df: int,
    hessian: list[list[float]],
    estimate_hill_coefficient: bool,
    warnings: list[str],
) -> tuple[dict[str, float | None], dict[str, float | None], dict[str, float | None]]:
    empty: dict[str, float | None] = {
        "standard_error": None,
        "ci_lower": None,
        "ci_upper": None,
        "ci_level": None,
    }
    if df < 1:
        warnings.append("Zero residual degrees of freedom; standard errors and confidence intervals were omitted.")
        return dict(empty), dict(empty), dict(empty)
    inv = invert_symmetric(hessian)
    if inv is None:
        warnings.append("Information matrix was singular; standard errors and confidence intervals were omitted.")
        return dict(empty), dict(empty), dict(empty)
    sigma2 = max(0.0, sse) / float(df)
    if not math.isfinite(sigma2):
        warnings.append("Residual variance was not finite; confidence intervals were omitted.")
        return dict(empty), dict(empty), dict(empty)
    t_crit = student_t_two_sided_quantile(df, _CI_LEVEL)

    se_emax = math.sqrt(sigma2 * inv[0][0])
    emax_unc: dict[str, float | None] = {
        "standard_error": se_emax,
        "ci_lower": emax - t_crit * se_emax,
        "ci_upper": emax + t_crit * se_emax,
        "ci_level": _CI_LEVEL,
    }
    se_eta = math.sqrt(sigma2 * inv[1][1])
    half_eta = min(t_crit * se_eta, 50.0)
    ec50_unc: dict[str, float | None] = {
        "standard_error": ec50 * se_eta,
        "ci_lower": ec50 * math.exp(-half_eta),
        "ci_upper": ec50 * math.exp(half_eta),
        "ci_level": _CI_LEVEL,
    }
    if not estimate_hill_coefficient:
        return emax_unc, ec50_unc, dict(empty)
    se_xi = math.sqrt(sigma2 * inv[2][2])
    half_xi = min(t_crit * se_xi, 50.0)
    hill_unc: dict[str, float | None] = {
        "standard_error": hill * se_xi,
        "ci_lower": hill * math.exp(-half_xi),
        "ci_upper": hill * math.exp(half_xi),
        "ci_level": _CI_LEVEL,
    }
    return emax_unc, ec50_unc, hill_unc


def _fit_linearized(
    *,
    fit_id: str,
    rows: tuple[DirectEffectPair, ...],
    concentrations: list[float],
    effects: list[float],
    concentration_unit: str,
    e0: float,
    e0_unit: str,
    direction: EffectDirection,
    effect_ids: tuple[str, ...],
    exposure_refs: tuple[str, ...],
) -> DirectEmaxFit:
    scale = max(1.0, abs(e0), max(abs(value) for value in effects))
    if any(abs(value - e0) <= 1e-12 * scale for value in effects):
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "An effect equals the fixed baseline, so the linear Emax rearrangement is undefined.",
        )
    deltas = [value - e0 for value in effects]
    if any(delta > 0 for delta in deltas) and any(delta < 0 for delta in deltas):
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Effects move both above and below the fixed baseline, so one Emax is not identifiable.",
        )
    span = max(deltas) - min(deltas)
    if span <= 1e-12 * max(1.0, max(abs(delta) for delta in deltas)):
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Effects do not vary with concentration, so Emax and EC50 are not identifiable.",
        )
    if max(concentrations) - min(concentrations) <= 0:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Concentrations do not vary, so EC50 is not identifiable.",
        )

    slope_intercept = _line([1.0 / value for value in concentrations], [1.0 / delta for delta in deltas])
    if slope_intercept is None:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Concentrations do not vary, so EC50 is not identifiable.",
        )
    slope, intercept = slope_intercept
    if not math.isfinite(slope) or not math.isfinite(intercept) or intercept == 0:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "Emax is not identifiable.",
        )
    emax = 1.0 / intercept
    ec50 = slope / intercept
    if not math.isfinite(emax) or not math.isfinite(ec50) or ec50 <= 0:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "EC50 is not positive, so the direct Emax fit is rejected.",
        )
    if direction is EffectDirection.INCREASE and emax <= 0:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "The declared direction is increase, but Emax is not positive. EC50 was not reported.",
        )
    if direction is EffectDirection.DECREASE and emax >= 0:
        return _failed(
            fit_id,
            DIRECT_EMAX,
            _METHOD_LINEARIZED,
            _LINEARIZED_LIMITATIONS,
            effect_ids,
            exposure_refs,
            rows,
            e0,
            e0_unit,
            "The declared direction is decrease, but Emax is not negative. EC50 was not reported.",
        )

    predicted_values = [
        emax_effect(concentration=value, e0=e0, emax=emax, ec50=ec50, hill=1.0) for value in concentrations
    ]
    residual = max(abs(pred - observed) for pred, observed in zip(predicted_values, effects, strict=True))
    status = DiagnosticStatus.PASS if residual <= 1e-6 * scale else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["An effect residual exceeds the noise-free tolerance."]
    estimates = (
        EffectParameter(
            fit_id=fit_id,
            name="E0",
            value=e0,
            unit=e0_unit,
            fixed=True,
            estimation_method="fixed",
        ),
        EffectParameter(
            fit_id=fit_id,
            name="Emax",
            value=emax,
            unit=e0_unit,
            fixed=False,
            estimation_method=_METHOD_LINEARIZED,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="EC50",
            value=ec50,
            unit=concentration_unit,
            fixed=False,
            estimation_method=_METHOD_LINEARIZED,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="hill",
            value=1.0,
            unit="1",
            fixed=True,
            estimation_method="fixed",
        ),
    )
    predicted = tuple(
        PredictedEffect(
            effect_id=row.effect_id,
            exposure_ref=row.exposure_ref,
            value=value,
            unit=e0_unit,
            derived_from=(f"fit:{fit_id}", f"exposure:{row.exposure_ref}"),
        )
        for row, value in zip(rows, predicted_values, strict=True)
    )
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(_LINEARIZED_LIMITATIONS),
        warnings=warnings,
        metrics=[
            Metric(name="n_points", value=float(len(rows))),
            Metric(name="max_abs_residual", value=residual, unit=e0_unit),
            Metric(name="Emax", value=emax, unit=e0_unit),
            Metric(name="EC50", value=ec50, unit=concentration_unit),
        ],
    )
    fit = _model_fit(
        fit_id,
        DIRECT_EMAX,
        _METHOD_LINEARIZED,
        status,
        converged=True,
        limitations=list(_LINEARIZED_LIMITATIONS),
    )
    return DirectEmaxFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        predicted=predicted,
        effect_ids=effect_ids,
        exposure_refs=exposure_refs,
        reason="",
        pairs=rows,
        e0=e0,
        e0_unit=e0_unit,
    )


def _logistic(z: float) -> float:
    if z >= 40.0:
        return 1.0
    if z <= -40.0:
        return 0.0
    if z >= 0.0:
        ez = math.exp(-z)
        return 1.0 / (1.0 + ez)
    ez = math.exp(z)
    return ez / (1.0 + ez)


def _line(xs: list[float], ys: list[float]) -> tuple[float, float] | None:
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0:
        return None
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    return slope, ybar - slope * xbar


def _failed(
    fit_id: str,
    spec: ModelSpec,
    method_tag: str,
    base_limitations: tuple[str, ...],
    effect_ids: tuple[str, ...],
    exposure_refs: tuple[str, ...],
    rows: tuple[DirectEffectPair, ...],
    e0: float,
    e0_unit: str,
    reason: str,
) -> DirectEmaxFit:
    stated = [*base_limitations, reason]
    warning = (
        "Emax, EC50, and hill were not estimated."
        if spec.model_id == SIGMOID_EMAX.model_id
        else "Emax and EC50 were not estimated."
    )
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=[warning],
    )
    fit = _model_fit(fit_id, spec, method_tag, DiagnosticStatus.FAIL, converged=False, limitations=stated)
    return DirectEmaxFit(
        fit=fit,
        estimates=(),
        diagnostic=diagnostic,
        predicted=(),
        effect_ids=effect_ids,
        exposure_refs=exposure_refs,
        reason=reason,
        pairs=rows,
        e0=e0,
        e0_unit=e0_unit,
    )


def _model_fit(
    fit_id: str,
    spec: ModelSpec,
    method_tag: str,
    status: DiagnosticStatus,
    *,
    converged: bool,
    limitations: list[str],
) -> ModelFit:
    return ModelFit(
        fit_id=fit_id,
        model_id=spec.model_id,
        dataset_id=fit_id,
        software_version=__version__,
        estimation_method=method_tag,
        context_of_use=spec.context_of_use,
        diagnostic_status=status,
        converged=converged,
        limitations=limitations,
    )
