"""EphysWeave endpoint attachment.

Reading an endpoint requires subject, session, timebase, and sync error.
This module does not shift clocks or interpolate onto a dose time.
"""

from pathlib import Path

from pharmacoweave.core.evidence import EvidenceRecord
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.integrations.guards import EphysEndpointRef, attach_evidence_file, read_ephys_endpoint


def read_neural_endpoint(record: EvidenceRecord) -> EphysEndpointRef:
    """Return an endpoint reference. The dose clock is not modified."""

    return read_ephys_endpoint(record)


def attach_neural_endpoints(study: PharmacologyStudy, path: str | Path) -> PharmacologyStudy:
    """Attach observed neural or behavioral endpoints that carry sync metadata."""

    return attach_evidence_file(study, path)
