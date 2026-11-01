"""PBPK is not calculated in this package.

PK-Sim is copyleft and is not imported. This package does not implement a second PBPK model.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError

STATUS = status("pbpk")


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process PBPK model."""

    raise PlannedCapabilityError(
        module="pbpk",
        phase="P4",
        detail=(
            "In-process PBPK is not implemented. PK-Sim is copyleft and is not imported or launched. "
            "This package does not reimplement PK-Sim."
        ),
    )
