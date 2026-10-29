"""Validation inventory and the research-use analysis report."""

from pharmacoweave.core.roadmap import status
from pharmacoweave.reporting.analysis_html import render_analysis_report
from pharmacoweave.reporting.html import render_validation_report

STATUS = status("reporting")

__all__ = ["STATUS", "render_analysis_report", "render_validation_report"]
