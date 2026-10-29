"""Diagnostic records and IV bolus figures.

A passing diagnostic still has to state limitations. Figures do not by
themselves mark a fit as acceptable.
"""

from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus
from pharmacoweave.core.roadmap import status

STATUS = status("diagnostics")

__all__ = ["DiagnosticRecord", "DiagnosticStatus", "STATUS"]
