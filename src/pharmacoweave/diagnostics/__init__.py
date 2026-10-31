"""Diagnostic records and the prior-versus-estimate prediction audit.

A passing record must state limitations. Goodness of fit alone cannot mark
a result as acceptable. Observed-versus-predicted and residual-versus-time
SVG figures are written by the result bundle. ``prediction_audit`` compares
predicted priors (TheraWeave, CADD, exploratory IVIVE) with estimated in vivo
parameters using fold-error metrics, without editing either side.
"""

from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus
from pharmacoweave.core.roadmap import status
from pharmacoweave.diagnostics.prediction_audit import (
    PREDICTION_AUDIT_LIMITATION,
    BiasDirection,
    EstimateContext,
    FoldMetrics,
    PredictionAuditEntry,
    PredictionAuditReport,
    PredictionAuditSummary,
    SkippedPredictionAudit,
    audit_evidence_predictions,
    audit_study_predictions,
    build_estimate_contexts,
    compute_fold_metrics,
    summarize_audit_entries,
)

STATUS = status("diagnostics")

__all__ = [
    "PREDICTION_AUDIT_LIMITATION",
    "STATUS",
    "BiasDirection",
    "DiagnosticRecord",
    "DiagnosticStatus",
    "EstimateContext",
    "FoldMetrics",
    "PredictionAuditEntry",
    "PredictionAuditReport",
    "PredictionAuditSummary",
    "SkippedPredictionAudit",
    "audit_evidence_predictions",
    "audit_study_predictions",
    "build_estimate_contexts",
    "compute_fold_metrics",
    "summarize_audit_entries",
]
