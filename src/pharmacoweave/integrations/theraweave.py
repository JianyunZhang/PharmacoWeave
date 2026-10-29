"""TheraWeave candidate attachment.

Predicted ADMET is stored as predicted evidence. It is not written into the
concentration table and it does not create a clearance estimate.
"""

from pathlib import Path

from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.integrations.guards import attach_evidence_file


def attach_candidate_bundle(study: PharmacologyStudy, path: str | Path) -> PharmacologyStudy:
    """Attach a TheraWeave evidence exchange after the project-boundary checks."""

    return attach_evidence_file(study, path)
