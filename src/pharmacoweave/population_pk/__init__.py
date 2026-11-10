"""Population PK is not calculated in this package.

nlmixr2 is copyleft and is not imported. The built-in executable backend is
only the MIT analytical reference engine. A caller-supplied executable is a
separate slot and does not add nlmixr2 to ``supported_backends()``.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError

STATUS = status("population_pk")


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process population model."""

    raise PlannedCapabilityError(
        module="population_pk",
        phase="P4",
        detail=(
            "In-process NLME is not implemented. nlmixr2 is copyleft and is not imported or launched. "
            "describe_engine('nlmixr2') records the license; supported_backends() does not include it."
        ),
    )
