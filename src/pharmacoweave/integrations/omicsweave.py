"""OmicsWeave evidence attachment.

Fluxes, transcripts, proteins, and pathway scores stay evidence records.
They are not drug concentrations or clearance estimates.
"""

from pathlib import Path

from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.integrations.guards import attach_evidence_file


def attach_biomarker_evidence(study: PharmacologyStudy, path: str | Path) -> PharmacologyStudy:
    """Attach an OmicsWeave evidence exchange without promoting it to PK data."""

    return attach_evidence_file(study, path)
