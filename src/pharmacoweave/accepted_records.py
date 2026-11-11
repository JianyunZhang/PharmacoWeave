"""Evaluate one already accepted record from a JSON request.

The request does not start a subprocess and does not interpret an engine-slot response.
"""

from collections.abc import Sequence
from typing import Any

from pydantic import BaseModel

from pharmacoweave.errors import AnalysisError
from pharmacoweave.pbpk.physiology_table import PhysiologyTable, validate_physiology_table
from pharmacoweave.pbpk.steady_state import steady_state_one_tissue
from pharmacoweave.pd_models.effect_compartment import simulate_effect_compartment
from pharmacoweave.pk_models.fit_three_compartment import recover_three_compartment_iv_bolus
from pharmacoweave.population_pk.results import load_population_result_table
from pharmacoweave.population_pk.two_stage import (
    IndividualParameterEstimate,
    summarize_log_parameters,
    summarize_weighted_log_parameters,
)
from pharmacoweave.qsp.inhibited_turnover import simulate_constant_inhibition_turnover
from pharmacoweave.qsp.turnover import simulate_linear_turnover

_FORBIDDEN = {"dosing_recommendation", "constitutes_dosing_recommendation", "engine_response"}


def evaluate_accepted_record(payload: dict[str, Any]) -> BaseModel | Sequence[BaseModel]:
    """Run one accepted function. Unknown kinds and dosing fields reject the request."""

    if _FORBIDDEN.intersection(payload) or payload.get("kind") == "engine_response":
        raise AnalysisError("An accepted record cannot carry a dosing recommendation or an engine response.")
    kind = payload.get("kind")
    body = {key: value for key, value in payload.items() if key != "kind"}
    if kind == "linear_turnover":
        return simulate_linear_turnover(**body)
    if kind == "inhibited_turnover":
        return simulate_constant_inhibition_turnover(**body)
    if kind == "effect_compartment":
        return simulate_effect_compartment(**body)
    if kind == "log_parameters":
        return summarize_log_parameters(_estimates(body))
    if kind == "weighted_log_parameters":
        return summarize_weighted_log_parameters(_estimates(body))
    if kind == "physiology_table":
        return validate_physiology_table(PhysiologyTable.model_validate(body))
    if kind == "steady_state_pbpk":
        table = PhysiologyTable.model_validate(body.pop("table"))
        return steady_state_one_tissue(table, **body)
    if kind == "population_result_table":
        rows = body.get("rows")
        if not isinstance(rows, list):
            raise AnalysisError("A population result table request needs a rows list.")
        return load_population_result_table(rows)
    if kind == "three_compartment_recovery":
        return recover_three_compartment_iv_bolus(**body)
    raise AnalysisError(f"Accepted record kind {kind!r} is not available.")


def _estimates(body: dict[str, Any]) -> list[IndividualParameterEstimate]:
    rows = body.get("estimates")
    if not isinstance(rows, list):
        raise AnalysisError("A parameter summary request needs an estimates list.")
    return [IndividualParameterEstimate.model_validate(row) for row in rows]
