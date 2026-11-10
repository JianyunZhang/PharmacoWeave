"""One linear turnover equation. Networks are not simulated.

``run_in_process`` still refuses a QSP network. ``simulate_linear_turnover``
evaluates one supplied closed form and does not estimate parameters.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError
from pharmacoweave.qsp.turnover import TurnoverSimulation, simulate_linear_turnover

STATUS = status("qsp")

__all__ = [
    "STATUS",
    "TurnoverSimulation",
    "run_in_process",
    "simulate_linear_turnover",
]


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process QSP model."""

    raise PlannedCapabilityError(
        module="qsp",
        phase="P4",
        detail=(
            "QSP networks are not implemented. simulate_linear_turnover evaluates one supplied "
            "linear turnover equation and does not estimate its parameters."
        ),
    )
