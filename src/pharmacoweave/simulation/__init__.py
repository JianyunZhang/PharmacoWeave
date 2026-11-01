"""Exploratory single-dose exposure simulation.

Importing this module does not run an analysis and does not start an external engine.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.simulation.exposure import (
    dosing_recommendation,
    human_starting_dose_from_simulation,
    refuse_unscoped_design,
    simulate_exposure,
)

STATUS = status("simulation")

__all__ = [
    "STATUS",
    "dosing_recommendation",
    "human_starting_dose_from_simulation",
    "refuse_unscoped_design",
    "simulate_exposure",
]
