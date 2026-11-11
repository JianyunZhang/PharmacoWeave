"""Export JSON Schema from the Pydantic contracts.

The Python models are the source of truth. ``schemas/`` is the checked-in
copy, and tests fail when the two diverge.
"""

import json
from pathlib import Path

from pydantic import BaseModel

from pharmacoweave.core.evidence import EvidenceExchange
from pharmacoweave.core.modeling import AnalysisRecords
from pharmacoweave.core.provenance import ResultBundleManifest
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.diagnostics.prediction_audit import PredictionAuditReport
from pharmacoweave.dmpk.records import InVitroAdmeSet
from pharmacoweave.ingestion.mapping import IngestMappingConfig
from pharmacoweave.pbpk.physiology_table import PhysiologyTableCheck
from pharmacoweave.pbpk.steady_state import SteadyStatePbpk
from pharmacoweave.pd_models.effect_compartment import EffectCompartmentSimulation
from pharmacoweave.pd_models.emax import PdAnalysisRecords
from pharmacoweave.pk_models.fit_three_compartment import ThreeCompartmentRecovery
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_nca.urine import UrineNcaConfig
from pharmacoweave.population_pk.results import PopulationResultTable
from pharmacoweave.population_pk.two_stage import LogParameterSummary
from pharmacoweave.qsp.concentration_turnover import ConcentrationTurnoverSimulation
from pharmacoweave.qsp.inhibited_turnover import InhibitedTurnoverSimulation
from pharmacoweave.qsp.turnover import TurnoverSimulation
from pharmacoweave.simulation.engine_protocol import CallerEngineRun, EngineResponse
from pharmacoweave.simulation.exposure import ExposureSimulationResult
from pharmacoweave.toxicology.records import ToxEndpointSet

CONTRACTS: dict[str, type[BaseModel]] = {
    "result_bundle.schema.json": ResultBundleManifest,
    "evidence_exchange.schema.json": EvidenceExchange,
    "pharmacology_study.schema.json": PharmacologyStudy,
    "analysis_records.schema.json": AnalysisRecords,
    "nca_config.schema.json": NcaConfig,
    "urine_nca_config.schema.json": UrineNcaConfig,
    "pd_records.schema.json": PdAnalysisRecords,
    "dmpk_records.schema.json": InVitroAdmeSet,
    "toxicology_records.schema.json": ToxEndpointSet,
    "prediction_audit.schema.json": PredictionAuditReport,
    "ingest_mapping.schema.json": IngestMappingConfig,
    "exposure_simulation.schema.json": ExposureSimulationResult,
    "engine_exchange.schema.json": EngineResponse,
    "caller_engine_run.schema.json": CallerEngineRun,
    "log_parameter_summary.schema.json": LogParameterSummary,
    "physiology_table_check.schema.json": PhysiologyTableCheck,
    "turnover_simulation.schema.json": TurnoverSimulation,
    "steady_state_pbpk.schema.json": SteadyStatePbpk,
    "population_result_table.schema.json": PopulationResultTable,
    "inhibited_turnover.schema.json": InhibitedTurnoverSimulation,
    "concentration_turnover.schema.json": ConcentrationTurnoverSimulation,
    "effect_compartment.schema.json": EffectCompartmentSimulation,
    "three_compartment_recovery.schema.json": ThreeCompartmentRecovery,
}


def schema_directory() -> Path:
    """Return the directory that holds the checked-in or packaged schemas."""

    packaged = Path(__file__).resolve().parent / "data" / "schemas"
    if (packaged / "result_bundle.schema.json").is_file():
        return packaged
    repo = Path(__file__).resolve().parents[2] / "schemas"
    if (repo / "result_bundle.schema.json").is_file():
        return repo
    raise FileNotFoundError(
        "PharmacoWeave JSON schemas were not found. Regenerate them with `python -m pharmacoweave.schema_export`."
    )


def export_schemas(directory: str | Path) -> None:
    """Write one JSON Schema file per public contract model."""

    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    for name, model in CONTRACTS.items():
        payload = json_schema(model)
        text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
        (destination / name).write_text(text, encoding="utf-8")


def json_schema(model: type[BaseModel]) -> dict[str, object]:
    """Return a JSON Schema that rejects unknown object fields."""

    schema = model.model_json_schema()
    _lock_objects(schema)
    return schema


def _lock_objects(node: object) -> None:
    if isinstance(node, dict):
        properties = node.get("properties")
        if node.get("type") == "object" and isinstance(properties, dict):
            node.setdefault("additionalProperties", False)
        for value in node.values():
            _lock_objects(value)
        return
    if isinstance(node, list):
        for value in node:
            _lock_objects(value)


def main() -> None:
    """Write schemas into the repository ``schemas`` directory."""

    repo_schemas = Path(__file__).resolve().parents[2] / "schemas"
    export_schemas(repo_schemas)


if __name__ == "__main__":
    main()
