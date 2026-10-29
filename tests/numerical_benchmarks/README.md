# Numerical benchmarks

These tests cover the first executable NCA and the one-compartment IV bolus fit. They are research-use checks. Passing them is not a bioanalytical validation and not a regulatory acceptance.

## NCA reference

`reference_linlog.py` is an independent transcription of:

- linear-up/log-down trapezoids, including a log trapezoid when the second concentration is lower and both are positive;
- IV bolus C0 by log-slope back-extrapolation of the first two declining quantified points when the first time is after zero;
- terminal λz by ordinary least squares of ln(concentration) on time, using points strictly after Tmax, at least 3 points, choosing the largest adjusted R², with values within `1e-6` treated as tied and the tie given to more points.

PKNCA and Phoenix are not installed and are not dependencies. The frozen digits in `test_nca_iv_bolus.py` come from that transcription on the 6-decimal quantified curve in `examples/data/iv_bolus_single_dose/`. The 24 h BLQ row is not an input to the reference.

## Tolerances

| Comparison | Quantity | Acceptance |
|---|---|---|
| Library vs frozen trapezoidal reference | AUC_last, AUC∞, fraction, Cmax, Tmax, λz, t½ | absolute `1e-9` |
| Library vs analytical target | AUC∞ vs Dose/CL = 5 | absolute `1e-4` |
| Library vs analytical target | λz vs 0.2 h⁻¹ | absolute `1e-5` |
| Library vs analytical target | t½ vs ln(2)/0.2 | absolute `1e-4` |
| NCA CL and Vz vs generating values | subject 1: 0.05 L/h, 0.25 L; subject 2: 0.06 L/h, 0.30 L | relative `1e-4` |
| Noise-free fit vs generating CL and V | same absolute parameters | relative `1e-4` |
| Proportional noise, σ = 0.02, seed 20261029 | CL and V | relative `0.10` |

The analytical window is wider than the trapezoidal window because the fixture rounds concentrations to 6 decimals. The calculator must match the trapezoids more tightly than Dose/CL. AUC∞ on this fixture is about 5.0000026, not exactly 5.

BLQ: the 24 h sample stays `BLQ` with an empty number. AUC_last must not match the value obtained by adding a linear trapezoid from 12 h down to zero at 24 h.

Reporting rules from `configs/nca_single_dose_iv.yaml`:

- fewer than 3 points after Tmax, or adjusted R² below 0.90: omit λz, t½, and AUC∞;
- extrapolated fraction above 0.20: keep λz and t½ when the regression itself passed, and omit AUC∞, CL, and Vz;
- duplicate times, an above-ULOQ sample, or an internal censored sample under `keep_as_censored`: omit the affected parameters and record the reason.

## Fit failures

A fit is not marked pass, and its estimates are omitted, when the profile has fewer than two quantified times, the optimizer does not converge, a parameter hits a bound, or the Jacobian is rank deficient. A parameter within a relative `1e-3` of an optimizer bound is treated as a bound hit, because the objective can flatten before `active_mask` is set.
