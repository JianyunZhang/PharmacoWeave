"""Validation reports.

The HTML writer inventories data and validation findings. It does not report
NCA parameters or model estimates.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.reporting.html import render_validation_report

STATUS = status("reporting")

__all__ = ["STATUS", "render_validation_report"]
