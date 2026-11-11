"""Population summaries of finished individual estimates, and a reader for caller-supplied tables.

nlmixr2 is copyleft and is not imported. ``run_in_process`` still refuses NLME.
``summarize_log_parameters`` does not refit data and does not weight by the standard error.
``summarize_weighted_log_parameters`` weights the log estimates by that standard error.
``load_population_result_table`` starts no subprocess.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.errors import PlannedCapabilityError
from pharmacoweave.population_pk.results import (
    PopulationResultRow,
    PopulationResultTable,
    load_population_result_table,
)
from pharmacoweave.population_pk.two_stage import (
    IndividualParameterEstimate,
    LogParameterSummary,
    summarize_log_parameters,
    summarize_weighted_log_parameters,
)

STATUS = status("population_pk")

__all__ = [
    "IndividualParameterEstimate",
    "LogParameterSummary",
    "PopulationResultRow",
    "PopulationResultTable",
    "STATUS",
    "load_population_result_table",
    "run_in_process",
    "summarize_log_parameters",
    "summarize_weighted_log_parameters",
]


def run_in_process(*_args: object, **_kwargs: object) -> None:
    """Refuse an in-process population model."""

    raise PlannedCapabilityError(
        module="population_pk",
        phase="P4",
        detail=(
            "In-process NLME is not implemented. nlmixr2 is copyleft and is not imported or launched. "
            "summarize_log_parameters and summarize_weighted_log_parameters summarize existing individual "
            "estimates and are not posteriors. "
            "describe_engine('nlmixr2') records the license; supported_backends() does not include it."
        ),
    )
