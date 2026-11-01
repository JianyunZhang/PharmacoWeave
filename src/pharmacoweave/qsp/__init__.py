"""Quantitative systems pharmacology models are not implemented.

Importing this module does not run an analysis.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError

STATUS = status("qsp")


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process QSP model."""

    raise PlannedCapabilityError(
        module="qsp",
        phase="P4",
        detail="QSP simulation is not implemented. The executable simulator covers analytical single-dose PK only.",
    )
