"""Validation reports and SVG figures.

The HTML writer inventories data and validation findings. NCA rows and IV bolus
fit figures are included only when that analysis was run.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.reporting.html import render_validation_report

STATUS = status("reporting")

__all__ = ["STATUS", "render_validation_report"]
