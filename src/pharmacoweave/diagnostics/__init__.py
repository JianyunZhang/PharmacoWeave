"""Diagnostic records.

A passing record must state limitations. Goodness of fit alone cannot mark
a result as acceptable. Observed-versus-predicted and residual-versus-time
SVG figures for an IV bolus fit are written by the result bundle. This
package does not estimate uncertainty intervals.
"""

from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus
from pharmacoweave.core.roadmap import status

STATUS = status("diagnostics")

__all__ = ["DiagnosticRecord", "DiagnosticStatus", "STATUS"]
