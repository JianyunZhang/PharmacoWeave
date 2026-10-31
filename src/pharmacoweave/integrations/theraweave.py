"""TheraWeave candidate attachment and prior-prediction audit.

Predicted ADMET is stored as predicted evidence. It is not written into the
concentration table and it does not create a clearance estimate.
``audit_candidate_predictions`` compares those predictions (and any
exploratory IVIVE) with PharmacoWeave in vivo estimates using fold-error
metrics; neither side is modified.
"""

from collections.abc import Sequence
from pathlib import Path

from pharmacoweave.core.evidence import EvidenceType
from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.diagnostics.prediction_audit import PredictionAuditReport, audit_study_predictions
from pharmacoweave.dmpk.physiology import IviveClearancePrediction
from pharmacoweave.integrations.guards import attach_evidence_file
from pharmacoweave.integrations.pharmacology_evidence import (
    PkFitResult,
    collect_run_evidence,
    evidence_from_ivive,
)


def attach_candidate_bundle(study: PharmacologyStudy, path: str | Path) -> PharmacologyStudy:
    """Attach a TheraWeave evidence exchange after the project-boundary checks."""

    return attach_evidence_file(study, path)


def audit_candidate_predictions(
    study: PharmacologyStudy,
    *,
    nca_parameters: Sequence[NcaParameterSet] | None = None,
    pk_fits: Sequence[PkFitResult] | None = None,
    ivive_predictions: Sequence[IviveClearancePrediction] | None = None,
    audit_id: str = "AUDIT-001",
) -> PredictionAuditReport:
    """Audit attached predicted priors and IVIVE predictions against in vivo estimates.

    Priors are the ``predicted`` records already attached to ``study`` plus
    one predicted record per IVIVE prediction. Posteriors are the
    ``estimated`` records built from ``nca_parameters`` and ``pk_fits``.
    """

    estimates = [
        record
        for record in collect_run_evidence(run_id=audit_id, nca_parameters=nca_parameters, pk_fits=pk_fits)
        if record.evidence_type is EvidenceType.ESTIMATED
    ]
    priors = [record for record in study.external_evidence if record.evidence_type is EvidenceType.PREDICTED]
    priors.extend(evidence_from_ivive(item, run_id=audit_id) for item in ivive_predictions or ())
    return audit_study_predictions(study, estimates, priors=priors, audit_id=audit_id)
