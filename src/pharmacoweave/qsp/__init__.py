"""Linear turnover equations. Networks are not simulated.

``run_in_process`` still refuses a QSP network. The turnover functions evaluate
supplied closed forms and do not estimate parameters.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError
from pharmacoweave.qsp.concentration_turnover import (
    ConcentrationTurnoverSimulation,
    simulate_concentration_driven_turnover,
)
from pharmacoweave.qsp.inhibited_turnover import (
    InhibitedTurnoverSimulation,
    simulate_constant_inhibition_turnover,
)
from pharmacoweave.qsp.turnover import TurnoverSimulation, simulate_linear_turnover

STATUS = status("qsp")

__all__ = [
    "STATUS",
    "ConcentrationTurnoverSimulation",
    "InhibitedTurnoverSimulation",
    "TurnoverSimulation",
    "run_in_process",
    "simulate_concentration_driven_turnover",
    "simulate_constant_inhibition_turnover",
    "simulate_linear_turnover",
]


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process QSP model."""

    raise PlannedCapabilityError(
        module="qsp",
        phase="P4",
        detail=(
            "QSP networks are not implemented. simulate_linear_turnover and "
            "simulate_constant_inhibition_turnover evaluate supplied closed forms and do not estimate parameters."
        ),
    )
