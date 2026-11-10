"""Two-stage summary of existing individual parameter estimates.

The summary is the geometric mean and the sample variance of the log parameters.
It does not refit the data and it does not produce a population posterior.
"""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Identifier, Record
from pharmacoweave.core.units import CLEARANCE, VOLUME, convert, parse_unit
from pharmacoweave.errors import AnalysisError

ParameterName = Literal["CL", "V", "CL/F", "V/F"]
_CANONICAL_UNIT = {"CL": "L/h", "CL/F": "L/h", "V": "L", "V/F": "L"}
_DIMENSION = {"CL": CLEARANCE, "CL/F": CLEARANCE, "V": VOLUME, "V/F": VOLUME}
_FAMILY = {"CL": "clearance", "CL/F": "clearance", "V": "volume", "V/F": "volume"}
_LIMITATIONS = (
    "This is a descriptive summary of individual estimates supplied by the caller.",
    "The original individual fit is unchanged. This function did not refit the data.",
    "The geometric mean is not an NLME posterior and is not a dosing recommendation.",
)


class IndividualParameterEstimate(Record):
    """One finished individual estimate. The standard error is stored and not used as a weight."""

    subject_id: Identifier
    parameter_name: ParameterName
    estimate: float
    unit: str = Field(min_length=1)
    standard_error: float | None = None


class LogParameterSummary(Record):
    """Geometric mean and log-scale sample variance for one parameter name."""

    parameter_name: ParameterName
    n_subjects: int = Field(ge=2)
    geometric_mean: float
    log_variance: float
    unit: str
    limitations: list[str]
    schema_version: Literal["0.1.0"] = "0.1.0"


def summarize_log_parameters(estimates: list[IndividualParameterEstimate]) -> list[LogParameterSummary]:
    """Summarize individual estimates on the log scale. A failed check returns no means."""

    if not estimates:
        raise AnalysisError("A two-stage summary requires individual estimates.")
    _reject_mixed_families(estimates)
    grouped: dict[str, list[IndividualParameterEstimate]] = {}
    for item in estimates:
        grouped.setdefault(item.parameter_name, []).append(item)
    summaries: list[LogParameterSummary] = []
    for name in sorted(grouped):
        rows = grouped[name]
        if len(rows) < 2:
            raise AnalysisError(f"{name} has fewer than two subjects, so no geometric mean is reported.")
        subjects = [item.subject_id for item in rows]
        if len(subjects) != len(set(subjects)):
            raise AnalysisError(f"{name} repeats a subject id, so no geometric mean is reported.")
        values = [_canonical_estimate(item) for item in rows]
        logs = [math.log(value) for value in values]
        mean_log = math.fsum(logs) / len(logs)
        variance = math.fsum((item - mean_log) ** 2 for item in logs) / (len(logs) - 1)
        summaries.append(
            LogParameterSummary(
                parameter_name=name,  # type: ignore[arg-type]
                n_subjects=len(rows),
                geometric_mean=math.exp(mean_log),
                log_variance=variance,
                unit=_CANONICAL_UNIT[name],
                limitations=list(_LIMITATIONS),
            )
        )
    return summaries


def _reject_mixed_families(estimates: list[IndividualParameterEstimate]) -> None:
    present = {item.parameter_name for item in estimates}
    for family in ("clearance", "volume"):
        names = sorted(name for name in present if _FAMILY[name] == family)
        if len(names) > 1:
            joined = " and ".join(names)
            raise AnalysisError(f"{joined} cannot be summarized together.")


def _canonical_estimate(item: IndividualParameterEstimate) -> float:
    if not math.isfinite(item.estimate) or item.estimate <= 0:
        raise AnalysisError(f"{item.subject_id} {item.parameter_name} must be a positive estimate.")
    dimensions = parse_unit(item.unit).dimensions
    if dimensions != _DIMENSION[item.parameter_name]:
        raise AnalysisError(f"{item.parameter_name} unit {item.unit!r} is not the parameter dimension.")
    return convert(item.estimate, item.unit, _CANONICAL_UNIT[item.parameter_name])
