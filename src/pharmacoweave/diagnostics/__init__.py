"""Diagnostic records.

A record can describe a future fit. This package does not draw
goodness-of-fit figures and cannot mark a fit as acceptable from a metric alone.
"""

from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus
from pharmacoweave.core.roadmap import status

STATUS = status("diagnostics")

__all__ = ["DiagnosticRecord", "DiagnosticStatus", "STATUS"]
