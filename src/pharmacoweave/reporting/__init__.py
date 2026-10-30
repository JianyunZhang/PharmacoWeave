"""Validation reports and SVG figures.

The HTML writer inventories data and validation findings. NCA rows, PK fit
tables/figures, and PD Emax curves are included only when that analysis was run.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.reporting.figures import (
    EffectCurveOverlay,
    EffectPlotPoint,
    PredictedPoint,
    SemilogOverlay,
    concentration_effect_svg,
    concentration_time_semilog_svg,
    concentration_time_svg,
    goodness_of_fit_svg,
    residual_diagnostics_svg,
    weighted_residual_svg,
)
from pharmacoweave.reporting.html import render_validation_report

STATUS = status("reporting")

__all__ = [
    "STATUS",
    "EffectCurveOverlay",
    "EffectPlotPoint",
    "PredictedPoint",
    "SemilogOverlay",
    "concentration_effect_svg",
    "concentration_time_semilog_svg",
    "concentration_time_svg",
    "goodness_of_fit_svg",
    "render_validation_report",
    "residual_diagnostics_svg",
    "weighted_residual_svg",
]
