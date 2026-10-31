"""Weighted nonlinear least squares for a two-compartment IV bolus.

Estimates ``CL`` (L/h), ``V1`` (L), ``Q`` (L/h), and ``V2`` (L) by curve
stripping (method of residuals) followed by multi-start Levenberg–Marquardt in
log-parameter space with exact analytical log-Jacobians.

Profiles with fewer than five quantified points, non-declining trajectories, or
monoexponential collapse (``alpha ≈ beta`` or unidentifiable peripheral phase)
return ``DiagnosticStatus.FAIL`` with no parameter estimates.
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
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint, iv_bolus_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import two_compartment_bolus_mg_l, two_compartment_macro_constants
from pharmacoweave.pk_models.fit_iv import ParameterUncertainty, screen_log_linear
from pharmacoweave.pk_models.specs import TWO_COMPARTMENT_IV_BOLUS
from pharmacoweave.version import __version__

_WEIGHTS = {
    "uniform": "weighted_nls:uniform",
    "1/y": "weighted_nls:1/y",
    "1/y^2": "weighted_nls:1/y^2",
}
_CI_LEVEL = 0.95
_MIN_ALPHA_OVER_BETA = 1.25
_MIN_FAST_FRACTION_AT_T1 = 0.03
_MIN_SLOW_FRACTION_AT_T0 = 0.02
_MIN_CURVATURE_RATIO = 1e-6
_LIMITATIONS = (
    "Assumes a single-dose two-compartment IV bolus with central elimination.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton log-parameter approximation.",
    "Research use only. The estimates are not a dosing recommendation.",
)


@dataclass(frozen=True)
class TwoCompartmentFit:
    """One subject-level 2-compartment IV bolus fit and its diagnostic record."""

    fit: ModelFit
    estimates: tuple[ParameterEstimate, ...]
    diagnostic: DiagnosticRecord
    observation_ids: tuple[str, ...]
    predicted_mg_l: tuple[float, ...]
    times_h: tuple[float, ...]
    observed_mg_l: tuple[float, ...]
    subject_id: str


def fit_two_compartment_iv_bolus(
    study: PharmacologyStudy,
    *,
    weight: str = "1/y^2",
) -> list[TwoCompartmentFit]:
    """Fit each single-dose IV bolus profile in ``study`` with a 2-compartment model."""

    return [
        fit_two_compartment_iv_bolus_profile(profile, weight=weight)
        for profile in iv_bolus_profiles(study)
    ]


def fit_two_compartment_iv_bolus_profile(
    profile: ConcentrationProfile,
    *,
    weight: str = "1/y^2",
) -> TwoCompartmentFit:
    """Fit one IV bolus profile; fail cleanly when biexponential kinetics are not supported."""

    method = _method(weight)
    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("fit_two_compartment_iv_bolus_profile requires route iv_bolus.")
    limitations = (*_LIMITATIONS, f"Weights are {weight} on quantified concentrations.")
    points = profile.points
    if len(points) < 5 or profile.dose_mg <= 0:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than five quantified points, or the dose is not positive.",
        )

    overall_screen = screen_log_linear(points)
    if (
        not overall_screen.identifiable
        or not overall_screen.declining
        or overall_screen.kel is None
        or overall_screen.intercept is None
    ):
        return _failed(profile, method=method, limitations=limitations, reason=overall_screen.reason)

    times = [p.time_h for p in points]
    observations = [p.concentration_mg_l for p in points]
    seeds = _initial_seeds(
        points,
        dose_mg=profile.dose_mg,
        overall_kel=overall_screen.kel,
        overall_intercept=overall_screen.intercept,
    )
    best: tuple[float, float, float, float, float, int, list[list[float]]] | None = None
    for cl0, v10, q0, v20 in seeds:
        candidate = _lm_two_compartment(
            times,
            observations,
            dose_mg=profile.dose_mg,
            cl=cl0,
            v1=v10,
            q=q0,
            v2=v20,
            weight=weight,
        )
        if candidate is None:
            continue
        if best is None or candidate[4] < best[4]:
            best = candidate

    if best is None:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason=(
                "The two-compartment weighted nonlinear least squares fit did not converge "
                "to an identifiable biexponential solution."
            ),
        )

    cl, v1, q, v2, weighted_sse, iterations, hessian = best
    a_coef, alpha, b_coef, beta = two_compartment_macro_constants(
        dose_mg=profile.dose_mg,
        cl_l_h=cl,
        v1_l=v1,
        q_l_h=q,
        v2_l=v2,
    )
    c0 = a_coef + b_coef
    c_first = a_coef * math.exp(-alpha * times[0]) + b_coef * math.exp(-beta * times[0])
    fast_frac_t1 = (a_coef * math.exp(-alpha * times[0])) / max(c_first, 1e-30)
    slow_frac_t0 = b_coef / max(c0, 1e-30)
    curvature = normalized_curvature_ratio(hessian)
    inv_hessian = invert_symmetric(hessian)
    if (
        beta <= 0
        or alpha / beta < _MIN_ALPHA_OVER_BETA
        or fast_frac_t1 < _MIN_FAST_FRACTION_AT_T1
        or slow_frac_t0 < _MIN_SLOW_FRACTION_AT_T0
        or curvature < _MIN_CURVATURE_RATIO
        or inv_hessian is None
    ):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason=(
                "Peripheral distribution is not separately identifiable from a monoexponential decline "
                f"(alpha/beta = {alpha / max(beta, 1e-30):.3f}, curvature ratio = {curvature:.2e})."
            ),
        )

    predicted = tuple(
        two_compartment_bolus_mg_l(dose_mg=profile.dose_mg, cl_l_h=cl, v1_l=v1, q_l_h=q, v2_l=v2, time_h=t_h)
        for t_h in times
    )
    max_relative = max(abs(pred - obs) / obs for pred, obs in zip(predicted, observations, strict=True))
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."]
    df = len(times) - 4
    uncertainties = _uncertainties_from_inverse(
        values=(cl, v1, q, v2),
        weighted_sse=weighted_sse,
        df=df,
        inv_hessian=inv_hessian,
        warnings=warnings,
    )
    metrics = [
        Metric(name="weighted_sse", value=weighted_sse),
        Metric(name="max_relative_residual", value=max_relative),
        Metric(name="alpha_over_beta_ratio", value=alpha / beta),
        Metric(name="curvature_ratio", value=curvature),
        Metric(name="n_points", value=float(len(times))),
        Metric(name="n_iterations", value=float(iterations)),
    ]
    return _succeeded(
        profile,
        method=method,
        limitations=limitations,
        status=status,
        cl=cl,
        v1=v1,
        q=q,
        v2=v2,
        uncertainties=uncertainties,
        predicted=predicted,
        warnings=warnings,
        metrics=metrics,
    )


def _initial_seeds(
    points: tuple[QuantifiedPoint, ...],
    *,
    dose_mg: float,
    overall_kel: float,
    overall_intercept: float,
) -> list[tuple[float, float, float, float]]:
    seeds: list[tuple[float, float, float, float]] = []
    n = len(points)
    tail_count = max(3, n // 3)
    tail_screen = screen_log_linear(points[-tail_count:])
    if (
        tail_screen.identifiable
        and tail_screen.declining
        and tail_screen.kel is not None
        and tail_screen.intercept is not None
    ):
        b_est = math.exp(tail_screen.intercept)
        beta_est = tail_screen.kel
        stripped: list[QuantifiedPoint] = []
        for idx, pt in enumerate(points[: n - tail_count]):
            diff = pt.concentration_mg_l - b_est * math.exp(-beta_est * pt.time_h)
            if diff > 1e-9:
                stripped.append(QuantifiedPoint(f"STRIP-{idx}", pt.time_h, diff))
        if len(stripped) >= 2:
            t1, y1 = stripped[0].time_h, stripped[0].concentration_mg_l
            t2, y2 = stripped[-1].time_h, stripped[-1].concentration_mg_l
            if t2 > t1 and y1 > y2 > 0:
                alpha_est = max(beta_est * 1.8, math.log(y1 / y2) / (t2 - t1))
                a_est = y1 * math.exp(alpha_est * t1)
                c0_est = a_est + b_est
                if c0_est > 0 and alpha_est > beta_est:
                    v1_s = dose_mg / c0_est
                    k21_s = (a_est * beta_est + b_est * alpha_est) / c0_est
                    if 0 < k21_s < alpha_est:
                        k10_s = (alpha_est * beta_est) / k21_s
                        k12_s = max(alpha_est + beta_est - k21_s - k10_s, 1e-4)
                        cl_s = k10_s * v1_s
                        q_s = k12_s * v1_s
                        v2_s = q_s / k21_s
                        if cl_s > 0 and v1_s > 0 and q_s > 0 and v2_s > 0:
                            seeds.append((cl_s, v1_s, q_s, v2_s))

    v_app = max(dose_mg / max(points[0].concentration_mg_l, math.exp(overall_intercept)), 1e-4)
    cl_app = max(overall_kel * v_app, 1e-5)
    for v1_mult in (0.5, 0.85, 1.1):
        v1 = v_app * v1_mult
        for cl_mult in (0.7, 1.2):
            cl = cl_app * cl_mult
            for q_ratio in (0.5, 2.5):
                q = cl * q_ratio
                for v2_ratio in (1.0, 3.0):
                    v2 = v1 * v2_ratio
                    seeds.append((cl, v1, q, v2))
    return seeds


def _macro_and_log_sensitivities(
    *,
    dose_mg: float,
    cl: float,
    v1: float,
    q: float,
    v2: float,
) -> tuple[float, float, float, float, list[tuple[float, float, float, float]]]:
    k10 = cl / v1
    k12 = q / v1
    k21 = q / v2
    s_rate = k10 + k12 + k21
    p_rate = k10 * k21
    disc_inner = (k10 + k12 - k21) * (k10 + k12 - k21) + 4.0 * k12 * k21
    r_span = math.sqrt(max(1e-30, disc_inner))
    alpha = 0.5 * (s_rate + r_span)
    beta = (2.0 * p_rate) / (s_rate + r_span)
    c0 = dose_mg / v1
    a_coef = c0 * (alpha - k21) / r_span
    b_coef = c0 - a_coef

    micro_log_dirs = (
        (k10, 0.0, 0.0, 0.0),          # d / d ln(CL)
        (-k10, -k12, 0.0, -c0),        # d / d ln(V1)
        (0.0, k12, k21, 0.0),          # d / d ln(Q)
        (0.0, 0.0, -k21, 0.0),         # d / d ln(V2)
    )
    sensitivities: list[tuple[float, float, float, float]] = []
    for dk10, dk12, dk21, dc0 in micro_log_dirs:
        ds = dk10 + dk12 + dk21
        dp = dk10 * k21 + k10 * dk21
        dr = (s_rate * ds - 2.0 * dp) / r_span
        d_alpha = 0.5 * (ds + dr)
        d_beta = 0.5 * (ds - dr)
        d_a = dc0 * (alpha - k21) / r_span + c0 * ((d_alpha - dk21) * r_span - (alpha - k21) * dr) / (r_span * r_span)
        d_b = dc0 - d_a
        sensitivities.append((d_a, d_alpha, d_b, d_beta))
    return a_coef, alpha, b_coef, beta, sensitivities


def _lm_two_compartment(
    times: list[float],
    observations: list[float],
    *,
    dose_mg: float,
    cl: float,
    v1: float,
    q: float,
    v2: float,
    weight: str,
) -> tuple[float, float, float, float, float, int, list[list[float]]] | None:
    def weight_of(obs: float) -> float:
        if weight == "uniform":
            return 1.0
        if weight == "1/y":
            return 1.0 / obs
        return 1.0 / (obs * obs)

    def sse(cl_v: float, v1_v: float, q_v: float, v2_v: float) -> float:
        a_c, alp, b_c, bet = two_compartment_macro_constants(
            dose_mg=dose_mg,
            cl_l_h=cl_v,
            v1_l=v1_v,
            q_l_h=q_v,
            v2_l=v2_v,
        )
        total = 0.0
        for t_h, obs in zip(times, observations, strict=True):
            if alp * t_h > 700.0 and bet * t_h > 700.0:
                pred = 0.0
            else:
                e_a = 0.0 if alp * t_h > 700.0 else math.exp(-alp * t_h)
                e_b = 0.0 if bet * t_h > 700.0 else math.exp(-bet * t_h)
                pred = a_c * e_a + b_c * e_b
            res = obs - pred
            total += weight_of(obs) * res * res
        return total

    def normal_equations(
        cl_v: float,
        v1_v: float,
        q_v: float,
        v2_v: float,
    ) -> tuple[list[list[float]], list[float]]:
        a_c, alp, b_c, bet, sens = _macro_and_log_sensitivities(
            dose_mg=dose_mg,
            cl=cl_v,
            v1=v1_v,
            q=q_v,
            v2=v2_v,
        )
        h = [[0.0, 0.0, 0.0, 0.0] for _ in range(4)]
        grad = [0.0, 0.0, 0.0, 0.0]
        for t_h, obs in zip(times, observations, strict=True):
            e_a = 0.0 if alp * t_h > 700.0 else math.exp(-alp * t_h)
            e_b = 0.0 if bet * t_h > 700.0 else math.exp(-bet * t_h)
            pred = a_c * e_a + b_c * e_b
            cols = [
                (da - a_c * t_h * d_alp) * e_a + (db - b_c * t_h * d_bet) * e_b
                for da, d_alp, db, d_bet in sens
            ]
            w = weight_of(obs)
            res = obs - pred
            for r in range(4):
                grad[r] += w * cols[r] * res
                for c in range(r, 4):
                    h[r][c] += w * cols[r] * cols[c]
        for r in range(4):
            for c in range(r):
                h[r][c] = h[c][r]
        return h, grad

    current = sse(cl, v1, q, v2)
    if not math.isfinite(current):
        return None
    lam = 1e-3
    converged = False
    iterations = 0
    for _step in range(90):
        iterations += 1
        hessian, grad = normal_equations(cl, v1, q, v2)
        if max(abs(g) for g in grad) <= 1e-9 * (1.0 + current):
            converged = True
            break
        accepted = False
        for _attempt in range(22):
            damped = [row[:] for row in hessian]
            for idx in range(4):
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
            if max_step > 2.0:
                factor = 2.0 / max_step
                step = [s * factor for s in step]
            u_cl = cl * math.exp(step[0])
            u_v1 = v1 * math.exp(step[1])
            u_q = q * math.exp(step[2])
            u_v2 = v2 * math.exp(step[3])
            if not all(1e-7 < val < 1e6 and math.isfinite(val) for val in (u_cl, u_v1, u_q, u_v2)):
                lam = min(lam * 8.0, 1e12)
                continue
            trial = sse(u_cl, u_v1, u_q, u_v2)
            if not math.isfinite(trial):
                lam = min(lam * 8.0, 1e12)
                continue
            if trial <= current * (1.0 + 1e-14):
                cl, v1, q, v2, current = u_cl, u_v1, u_q, u_v2, trial
                lam = max(lam / 3.0, 1e-14)
                accepted = True
                break
            lam = min(lam * 8.0, 1e12)
        if converged or not accepted:
            break

    final_h, _final_g = normal_equations(cl, v1, q, v2)
    if not converged or not (cl > 0 and v1 > 0 and q > 0 and v2 > 0 and math.isfinite(current)):
        return None
    return cl, v1, q, v2, current, iterations, final_h


def _uncertainties_from_inverse(
    *,
    values: tuple[float, float, float, float],
    weighted_sse: float,
    df: int,
    inv_hessian: list[list[float]],
    warnings: list[str],
) -> tuple[ParameterUncertainty, ParameterUncertainty, ParameterUncertainty, ParameterUncertainty]:
    empty = (
        ParameterUncertainty(),
        ParameterUncertainty(),
        ParameterUncertainty(),
        ParameterUncertainty(),
    )
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
    return out[0], out[1], out[2], out[3]


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
    cl: float,
    v1: float,
    q: float,
    v2: float,
    uncertainties: tuple[ParameterUncertainty, ParameterUncertainty, ParameterUncertainty, ParameterUncertainty],
    predicted: tuple[float, ...],
    warnings: list[str],
    metrics: list[Metric],
) -> TwoCompartmentFit:
    fit_id = _fit_id(profile)
    cl_unc, v1_unc, q_unc, v2_unc = uncertainties
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(limitations),
        warnings=warnings,
        metrics=metrics,
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=TWO_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=TWO_COMPARTMENT_IV_BOLUS.context_of_use,
        diagnostic_status=status,
        converged=True,
        limitations=list(limitations),
    )
    estimates = (
        ParameterEstimate(
            fit_id=fit_id,
            parameter="CL",
            semantic=ParameterSemantic.CL,
            value=cl,
            unit="L/h",
            estimation_method=method,
            route=Route.IV_BOLUS,
            standard_error=cl_unc.standard_error,
            ci_lower=cl_unc.ci_lower,
            ci_upper=cl_unc.ci_upper,
            ci_level=cl_unc.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V1",
            semantic=ParameterSemantic.V1,
            value=v1,
            unit="L",
            estimation_method=method,
            route=Route.IV_BOLUS,
            standard_error=v1_unc.standard_error,
            ci_lower=v1_unc.ci_lower,
            ci_upper=v1_unc.ci_upper,
            ci_level=v1_unc.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="Q",
            semantic=ParameterSemantic.Q,
            value=q,
            unit="L/h",
            estimation_method=method,
            route=Route.IV_BOLUS,
            standard_error=q_unc.standard_error,
            ci_lower=q_unc.ci_lower,
            ci_upper=q_unc.ci_upper,
            ci_level=q_unc.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V2",
            semantic=ParameterSemantic.V2,
            value=v2,
            unit="L",
            estimation_method=method,
            route=Route.IV_BOLUS,
            standard_error=v2_unc.standard_error,
            ci_lower=v2_unc.ci_lower,
            ci_upper=v2_unc.ci_upper,
            ci_level=v2_unc.ci_level,
        ),
    )
    return _pack(profile, fit, estimates, diagnostic, predicted)


def _failed(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    reason: str,
) -> TwoCompartmentFit:
    fit_id = _fit_id(profile)
    stated = [*limitations, reason]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=["CL, V1, Q, and V2 were not estimated."],
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=TWO_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=TWO_COMPARTMENT_IV_BOLUS.context_of_use,
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
) -> TwoCompartmentFit:
    return TwoCompartmentFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        observation_ids=tuple(point.observation_id for point in profile.points),
        predicted_mg_l=predicted,
        times_h=tuple(point.time_h for point in profile.points),
        observed_mg_l=tuple(point.concentration_mg_l for point in profile.points),
        subject_id=profile.subject_id,
    )


def _fit_id(profile: ConcentrationProfile) -> str:
    return f"fit-2c-iv-{profile.subject_id}-{profile.analyte_id}"


def _dataset_id(profile: ConcentrationProfile) -> str:
    return f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}"
