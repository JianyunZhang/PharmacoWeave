"""Weighted nonlinear least squares for single-dose 1-compartment extravascular absorption.

Estimates ``ka`` (1/h), apparent clearance ``CL/F`` (L/h), and apparent central
volume ``V/F`` (L) by multi-start Levenberg–Marquardt in log-parameter space
under the canonical assumption ``ka > ke = (CL/F)/(V/F)``.

Profiles with fewer than four quantified points, no rising absorption-phase
samples before ``Tmax``, fewer than two post-peak elimination samples, or
non-identifiable ``ka ≈ ke`` / flip-flop confluence return
``DiagnosticStatus.FAIL`` with no parameter estimates.
"""

import math
from dataclasses import dataclass

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import (
    DiagnosticRecord,
    DiagnosticStatus,
    Metric,
    ModelFit,
    ParameterEstimate,
    ParameterSemantic,
)
from pharmacoweave.core.numerics import (
    invert_symmetric,
    normalized_curvature_ratio,
    solve_linear_system,
    student_t_two_sided_quantile,
)
from pharmacoweave.core.profiles import ConcentrationProfile, extravascular_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import bateman_mg_l
from pharmacoweave.pk_models.fit_iv import ParameterUncertainty, screen_log_linear
from pharmacoweave.pk_models.specs import ONE_COMPARTMENT_EXTRAVASCULAR
from pharmacoweave.version import __version__

_EXTRAVASCULAR_ROUTES = {Route.ORAL, Route.SUBCUTANEOUS, Route.INTRAMUSCULAR}
_WEIGHTS = {
    "uniform": "weighted_nls:uniform",
    "1/y": "weighted_nls:1/y",
    "1/y^2": "weighted_nls:1/y^2",
}
_CI_LEVEL = 0.95
_MIN_KA_OVER_KE = 1.15
_MIN_CURVATURE_RATIO = 1e-6
_LIMITATIONS = (
    "Assumes single-dose one-compartment first-order absorption without lag time.",
    "Identifies the canonical ka > ke branch; flip-flop kinetics cannot be distinguished without an IV arm.",
    "Reports apparent CL/F and V/F because bioavailability F is unknown in a single-arm extravascular design.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton log-parameter approximation.",
    "Research use only. The estimates are not a dosing recommendation.",
)


@dataclass(frozen=True)
class ExtravascularFit:
    """One subject-level 1-compartment extravascular fit and its diagnostic record."""

    fit: ModelFit
    estimates: tuple[ParameterEstimate, ...]
    diagnostic: DiagnosticRecord
    observation_ids: tuple[str, ...]
    predicted_mg_l: tuple[float, ...]
    times_h: tuple[float, ...]
    observed_mg_l: tuple[float, ...]
    subject_id: str
    route: Route


def fit_one_compartment_extravascular(
    study: PharmacologyStudy,
    *,
    weight: str = "1/y^2",
) -> list[ExtravascularFit]:
    """Fit each single-dose extravascular profile in ``study`` by weighted NLS."""

    return [
        fit_one_compartment_extravascular_profile(profile, weight=weight)
        for profile in extravascular_profiles(study)
    ]


def fit_one_compartment_extravascular_profile(
    profile: ConcentrationProfile,
    *,
    weight: str = "1/y^2",
) -> ExtravascularFit:
    """Fit one single-dose extravascular profile; withhold estimates when unidentifiable."""

    method = _method(weight)
    if profile.route not in _EXTRAVASCULAR_ROUTES:
        raise AnalysisError(
            "fit_one_compartment_extravascular_profile requires route oral, subcutaneous, or intramuscular."
        )
    limitations = (*_LIMITATIONS, f"Weights are {weight} on quantified concentrations.")
    points = profile.points
    if len(points) < 4 or profile.dose_mg <= 0:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than four quantified points, or the dose is not positive.",
        )
    if any(p.time_h <= 0 or p.concentration_mg_l <= 0 for p in points):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="All quantified extravascular points must have positive post-dose time and concentration.",
        )

    cmax_idx = max(range(len(points)), key=lambda idx: (points[idx].concentration_mg_l, -points[idx].time_h))
    if cmax_idx < 1 or points[0].concentration_mg_l >= points[cmax_idx].concentration_mg_l * (1.0 - 1e-6):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="No rising absorption-phase points precede Tmax, so ka cannot be identified.",
        )
    post_peak = points[cmax_idx:]
    if len(points) - 1 - cmax_idx < 2:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than two post-peak elimination points follow Tmax.",
        )
    tail_points = points[max(cmax_idx + 1, len(points) - 3) :]
    if len(tail_points) < 3:
        tail_points = post_peak
    tail_screen = screen_log_linear(tail_points)
    if (
        not tail_screen.identifiable
        or not tail_screen.declining
        or tail_screen.kel is None
        or tail_screen.intercept is None
    ):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Post-peak concentrations do not define a declining elimination phase.",
        )

    times = [p.time_h for p in points]
    observations = [p.concentration_mg_l for p in points]
    seeds = _initial_seeds(
        dose_mg=profile.dose_mg,
        tmax_h=points[cmax_idx].time_h,
        cmax_mg_l=points[cmax_idx].concentration_mg_l,
        kel_tail=tail_screen.kel,
    )
    best: tuple[float, float, float, float, int, list[list[float]]] | None = None
    for ka0, clf0, vf0 in seeds:
        candidate = _lm_bateman(
            times,
            observations,
            dose_mg=profile.dose_mg,
            ka=ka0,
            clf=clf0,
            vf=vf0,
            weight=weight,
        )
        if candidate is None:
            continue
        if best is None or candidate[3] < best[3]:
            best = candidate

    if best is None:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason=(
                "The weighted nonlinear least squares fit did not converge to an identifiable "
                "ka > ke solution (possible ka ≈ ke confluence or flip-flop non-identifiability)."
            ),
        )

    ka, clf, vf, weighted_sse, iterations, hessian = best
    kel = clf / vf
    curvature = normalized_curvature_ratio(hessian)
    inv_hessian = invert_symmetric(hessian)
    if (
        ka / kel < _MIN_KA_OVER_KE
        or ka * times[0] > 18.0
        or curvature < _MIN_CURVATURE_RATIO
        or inv_hessian is None
    ):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason=(
                "Absorption and elimination rates are not separately identifiable "
                f"(ka/ke = {ka / kel:.3f}, curvature ratio = {curvature:.2e})."
            ),
        )

    predicted = tuple(
        bateman_mg_l(dose_mg=profile.dose_mg, ka_h=ka, clf_l_h=clf, vf_l=vf, time_h=t_h) for t_h in times
    )
    max_relative = max(abs(pred - obs) / obs for pred, obs in zip(predicted, observations, strict=True))
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."]
    df = len(times) - 3
    uncertainties = _uncertainties_from_inverse(
        values=(ka, clf, vf),
        weighted_sse=weighted_sse,
        df=df,
        inv_hessian=inv_hessian,
        warnings=warnings,
    )
    metrics = [
        Metric(name="weighted_sse", value=weighted_sse),
        Metric(name="max_relative_residual", value=max_relative),
        Metric(name="ka_over_ke_ratio", value=ka / kel),
        Metric(name="curvature_ratio", value=curvature),
        Metric(name="n_points", value=float(len(times))),
        Metric(name="n_iterations", value=float(iterations)),
    ]
    return _succeeded(
        profile,
        method=method,
        limitations=limitations,
        status=status,
        ka=ka,
        clf=clf,
        vf=vf,
        uncertainties=uncertainties,
        predicted=predicted,
        warnings=warnings,
        metrics=metrics,
    )


def _initial_seeds(
    *,
    dose_mg: float,
    tmax_h: float,
    cmax_mg_l: float,
    kel_tail: float,
) -> list[tuple[float, float, float]]:
    kel0 = max(kel_tail, 1e-3)
    base_vf = max(dose_mg / max(cmax_mg_l * 1.5, 1e-9), 1e-4)
    seeds: list[tuple[float, float, float]] = []
    for kel_mult in (0.85, 1.0, 1.2):
        kel = kel0 * kel_mult
        approx_ka_from_tmax = max(kel * 1.6, 1.5 / max(tmax_h, 0.05))
        for ka_mult in (0.6, 1.0, 2.0, 4.0):
            ka = max(kel * 1.4, approx_ka_from_tmax * ka_mult)
            vf = base_vf
            clf = kel * vf
            c_at_tmax = bateman_mg_l(dose_mg=dose_mg, ka_h=ka, clf_l_h=clf, vf_l=vf, time_h=tmax_h)
            if c_at_tmax > 1e-12 and cmax_mg_l > 0:
                scale = c_at_tmax / cmax_mg_l
                vf = max(vf * scale, 1e-6)
                clf = kel * vf
            seeds.append((ka, clf, vf))
    return seeds


def _bateman_and_log_jacobian(
    *,
    dose_mg: float,
    ka: float,
    clf: float,
    vf: float,
    time_h: float,
) -> tuple[float, float, float, float]:
    kel = clf / vf
    x = (ka - kel) * time_h
    exp_kel = math.exp(-kel * time_h)
    if abs(x) <= 1e-5:
        x2 = x * x
        x3 = x2 * x
        g_val = 1.0 - 0.5 * x + x2 / 6.0 - x3 / 24.0
        dg_val = -0.5 + x / 3.0 - x2 / 8.0 + x3 / 30.0
    else:
        g_val = -math.expm1(-x) / x
        dg_val = (math.exp(-x) * (1.0 + x) - 1.0) / (x * x)
    prefactor = (dose_mg / vf) * (ka * time_h) * exp_kel
    conc = prefactor * g_val
    col_ka = conc + prefactor * dg_val * (ka * time_h)
    d_kel = -time_h * (conc + prefactor * dg_val)
    col_clf = kel * d_kel
    col_vf = -conc - col_clf
    return conc, col_ka, col_clf, col_vf


def _lm_bateman(
    times: list[float],
    observations: list[float],
    *,
    dose_mg: float,
    ka: float,
    clf: float,
    vf: float,
    weight: str,
) -> tuple[float, float, float, float, int, list[list[float]]] | None:
    def weight_of(obs: float) -> float:
        if weight == "uniform":
            return 1.0
        if weight == "1/y":
            return 1.0 / obs
        return 1.0 / (obs * obs)

    def sse(ka_v: float, clf_v: float, vf_v: float) -> float:
        total = 0.0
        for t_h, obs in zip(times, observations, strict=True):
            pred = bateman_mg_l(dose_mg=dose_mg, ka_h=ka_v, clf_l_h=clf_v, vf_l=vf_v, time_h=t_h)
            res = obs - pred
            total += weight_of(obs) * res * res
        return total

    def normal_equations(ka_v: float, clf_v: float, vf_v: float) -> tuple[list[list[float]], list[float]]:
        h = [[0.0, 0.0, 0.0] for _ in range(3)]
        grad = [0.0, 0.0, 0.0]
        for t_h, obs in zip(times, observations, strict=True):
            pred, j0, j1, j2 = _bateman_and_log_jacobian(
                dose_mg=dose_mg,
                ka=ka_v,
                clf=clf_v,
                vf=vf_v,
                time_h=t_h,
            )
            cols = (j0, j1, j2)
            w = weight_of(obs)
            res = obs - pred
            for r in range(3):
                grad[r] += w * cols[r] * res
                for c in range(r, 3):
                    h[r][c] += w * cols[r] * cols[c]
        for r in range(3):
            for c in range(r):
                h[r][c] = h[c][r]
        return h, grad

    current = sse(ka, clf, vf)
    if not math.isfinite(current):
        return None
    lam = 1e-3
    converged = False
    iterations = 0
    for _step in range(80):
        iterations += 1
        hessian, grad = normal_equations(ka, clf, vf)
        if max(abs(g) for g in grad) <= 1e-9 * (1.0 + current):
            converged = True
            break
        accepted = False
        for _attempt in range(22):
            damped = [row[:] for row in hessian]
            for idx in range(3):
                scale = damped[idx][idx] if damped[idx][idx] > 1e-18 else 1.0
                damped[idx][idx] += lam * scale
            step = solve_linear_system(damped, grad)
            if step is None:
                lam = min(lam * 8.0, 1e12)
                continue
            max_step = max(abs(s) for s in step)
            if max_step <= 1e-10:
                converged = True
                accepted = True
                break
            if max_step > 2.5:
                factor = 2.5 / max_step
                step = [s * factor for s in step]
            u_ka = ka * math.exp(step[0])
            u_clf = clf * math.exp(step[1])
            u_vf = vf * math.exp(step[2])
            if not (math.isfinite(u_ka) and math.isfinite(u_clf) and math.isfinite(u_vf)):
                lam = min(lam * 8.0, 1e12)
                continue
            # Keep solver on the canonical ka >= ke branch to avoid flip-flop jumps
            if u_ka < (u_clf / u_vf) * 1.0001 or u_ka > 1e4 or (u_clf / u_vf) < 1e-6:
                lam = min(lam * 8.0, 1e12)
                continue
            trial = sse(u_ka, u_clf, u_vf)
            if not math.isfinite(trial):
                lam = min(lam * 8.0, 1e12)
                continue
            if trial <= current * (1.0 + 1e-14):
                ka, clf, vf, current = u_ka, u_clf, u_vf, trial
                lam = max(lam / 3.0, 1e-14)
                accepted = True
                break
            lam = min(lam * 8.0, 1e12)
        if converged or not accepted:
            break

    final_h, _final_g = normal_equations(ka, clf, vf)
    if not converged or not (ka > 0 and clf > 0 and vf > 0 and math.isfinite(current)):
        return None
    return ka, clf, vf, current, iterations, final_h


def _uncertainties_from_inverse(
    *,
    values: tuple[float, float, float],
    weighted_sse: float,
    df: int,
    inv_hessian: list[list[float]],
    warnings: list[str],
) -> tuple[ParameterUncertainty, ParameterUncertainty, ParameterUncertainty]:
    empty = (ParameterUncertainty(), ParameterUncertainty(), ParameterUncertainty())
    if df < 1:
        warnings.append("Residual degrees of freedom did not support confidence intervals.")
        return empty
    sigma2 = max(0.0, weighted_sse) / float(df)
    if not math.isfinite(sigma2):
        warnings.append("Residual variance was not finite; confidence intervals were omitted.")
        return empty
    t_crit = student_t_two_sided_quantile(df, _CI_LEVEL)
    out: list[ParameterUncertainty] = []
    for idx, val in enumerate(values):
        se_log = math.sqrt(max(0.0, sigma2 * inv_hessian[idx][idx]))
        half_width = min(t_crit * se_log, 50.0)
        out.append(
            ParameterUncertainty(
                standard_error=val * se_log,
                ci_lower=val * math.exp(-half_width),
                ci_upper=val * math.exp(half_width),
                ci_level=_CI_LEVEL,
            )
        )
    return out[0], out[1], out[2]


def _method(weight: str) -> str:
    try:
        return _WEIGHTS[weight]
    except KeyError as exc:
        allowed = ", ".join(sorted(_WEIGHTS))
        raise AnalysisError(f"weight must be one of {allowed}.") from exc


def _succeeded(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    status: DiagnosticStatus,
    ka: float,
    clf: float,
    vf: float,
    uncertainties: tuple[ParameterUncertainty, ParameterUncertainty, ParameterUncertainty],
    predicted: tuple[float, ...],
    warnings: list[str],
    metrics: list[Metric],
) -> ExtravascularFit:
    fit_id = _fit_id(profile)
    ka_unc, clf_unc, vf_unc = uncertainties
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(limitations),
        warnings=warnings,
        metrics=metrics,
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_EXTRAVASCULAR.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_EXTRAVASCULAR.context_of_use,
        diagnostic_status=status,
        converged=True,
        limitations=list(limitations),
    )
    estimates = (
        ParameterEstimate(
            fit_id=fit_id,
            parameter="ka",
            semantic=ParameterSemantic.KA,
            value=ka,
            unit="1/h",
            estimation_method=method,
            route=profile.route,
            standard_error=ka_unc.standard_error,
            ci_lower=ka_unc.ci_lower,
            ci_upper=ka_unc.ci_upper,
            ci_level=ka_unc.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="CL/F",
            semantic=ParameterSemantic.CL_OVER_F,
            value=clf,
            unit="L/h",
            estimation_method=method,
            route=profile.route,
            standard_error=clf_unc.standard_error,
            ci_lower=clf_unc.ci_lower,
            ci_upper=clf_unc.ci_upper,
            ci_level=clf_unc.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V/F",
            semantic=ParameterSemantic.V_OVER_F,
            value=vf,
            unit="L",
            estimation_method=method,
            route=profile.route,
            standard_error=vf_unc.standard_error,
            ci_lower=vf_unc.ci_lower,
            ci_upper=vf_unc.ci_upper,
            ci_level=vf_unc.ci_level,
        ),
    )
    return _pack(profile, fit, estimates, diagnostic, predicted)


def _failed(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    reason: str,
) -> ExtravascularFit:
    fit_id = _fit_id(profile)
    stated = [*limitations, reason]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=["ka, CL/F, and V/F were not estimated."],
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_EXTRAVASCULAR.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_EXTRAVASCULAR.context_of_use,
        diagnostic_status=DiagnosticStatus.FAIL,
        converged=False,
        limitations=stated,
    )
    return _pack(profile, fit, (), diagnostic, ())


def _pack(
    profile: ConcentrationProfile,
    fit: ModelFit,
    estimates: tuple[ParameterEstimate, ...],
    diagnostic: DiagnosticRecord,
    predicted: tuple[float, ...],
) -> ExtravascularFit:
    return ExtravascularFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        observation_ids=tuple(point.observation_id for point in profile.points),
        predicted_mg_l=predicted,
        times_h=tuple(point.time_h for point in profile.points),
        observed_mg_l=tuple(point.concentration_mg_l for point in profile.points),
        subject_id=profile.subject_id,
        route=profile.route,
    )


def _fit_id(profile: ConcentrationProfile) -> str:
    return f"fit-1c-ev-{profile.subject_id}-{profile.analyte_id}"


def _dataset_id(profile: ConcentrationProfile) -> str:
    return f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}"
