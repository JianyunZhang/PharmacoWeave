"""Explicit physiology checks and one algebraic liver-plus-tissue steady state.

PK-Sim is copyleft and is not imported. ``run_in_process`` still refuses a
dynamic whole-body PBPK model.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError
from pharmacoweave.pbpk.physiology_table import (
    PhysiologyTable,
    PhysiologyTableCheck,
    validate_physiology_table,
)
from pharmacoweave.pbpk.steady_state import SteadyStatePbpk, steady_state_one_tissue

STATUS = status("pbpk")

__all__ = [
    "PhysiologyTable",
    "PhysiologyTableCheck",
    "STATUS",
    "SteadyStatePbpk",
    "run_in_process",
    "steady_state_one_tissue",
    "validate_physiology_table",
]


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process PBPK model."""

    raise PlannedCapabilityError(
        module="pbpk",
        phase="P4",
        detail=(
            "Dynamic whole-body PBPK is not implemented. PK-Sim is copyleft and is not imported or launched. "
            "validate_physiology_table and steady_state_one_tissue use caller-supplied numbers only."
        ),
    )
