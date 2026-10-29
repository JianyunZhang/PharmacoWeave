"""External PK engine boundary.

NONMEM, nlmixr2, PK-Sim, and MoBi are not imported. nlmixr2 and parts of the
Open Systems Pharmacology stack use copyleft licenses, so a future adapter
has to exchange files through a separate process and record the engine
version in the result bundle.
"""

from pharmacoweave.errors import PlannedCapabilityError


def supported_backends() -> list[str]:
    """Return engine names this package can execute. The list is empty."""

    return []


def run_backend(name: str) -> None:
    """Refuse to launch an external engine from this version."""

    raise PlannedCapabilityError(
        module="integrations.external_backends",
        phase="P2",
        detail=(
            f"No adapter is registered for {name!r}. Supported backends: none. "
            "An adapter must run out of process and record the engine version and license."
        ),
    )
