"""Checked-in JSON Schemas stay aligned with the Pydantic models."""

import json

import jsonschema

from pharmacoweave.core.provenance import ResultBundleManifest
from pharmacoweave.diagnostics.prediction_audit import PredictionAuditReport
from pharmacoweave.dmpk.records import InVitroAdmeSet
from pharmacoweave.ingestion.mapping import IngestMappingConfig
from pharmacoweave.pbpk.physiology_table import PhysiologyTableCheck
from pharmacoweave.pbpk.steady_state import SteadyStatePbpk
from pharmacoweave.pd_models.effect_compartment import EffectCompartmentSimulation
from pharmacoweave.pd_models.emax import PdAnalysisRecords
from pharmacoweave.pk_models.fit_three_compartment import ThreeCompartmentRecovery
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.population_pk.results import PopulationResultTable
from pharmacoweave.population_pk.two_stage import LogParameterSummary
from pharmacoweave.qsp.inhibited_turnover import InhibitedTurnoverSimulation
from pharmacoweave.qsp.turnover import TurnoverSimulation
from pharmacoweave.schema_export import CONTRACTS, json_schema
from pharmacoweave.simulation.engine_protocol import EngineResponse
from pharmacoweave.simulation.exposure import ExposureSimulationResult
from pharmacoweave.toxicology.records import ToxEndpointSet
from pharmacoweave.version import SCHEMA_VERSION


def test_schema_version_defaults_match_the_package_constant() -> None:
    assert ResultBundleManifest.model_fields["schema_version"].default == SCHEMA_VERSION
    assert NcaConfig.model_fields["schema_version"].default == SCHEMA_VERSION
    assert PdAnalysisRecords.model_fields["schema_version"].default == SCHEMA_VERSION
    assert InVitroAdmeSet.model_fields["schema_version"].default == SCHEMA_VERSION
    assert ToxEndpointSet.model_fields["schema_version"].default == SCHEMA_VERSION
    assert PredictionAuditReport.model_fields["schema_version"].default == SCHEMA_VERSION
    assert IngestMappingConfig.model_fields["schema_version"].default == SCHEMA_VERSION
    assert ExposureSimulationResult.model_fields["schema_version"].default == "0.2.2"
    assert EngineResponse.model_fields["schema_version"].default == "0.2.2"
    assert LogParameterSummary.model_fields["schema_version"].default == SCHEMA_VERSION
    assert PhysiologyTableCheck.model_fields["schema_version"].default == SCHEMA_VERSION
    assert TurnoverSimulation.model_fields["schema_version"].default == SCHEMA_VERSION
    assert SteadyStatePbpk.model_fields["schema_version"].default == SCHEMA_VERSION
    assert PopulationResultTable.model_fields["schema_version"].default == SCHEMA_VERSION
    assert InhibitedTurnoverSimulation.model_fields["schema_version"].default == SCHEMA_VERSION
    assert EffectCompartmentSimulation.model_fields["schema_version"].default == SCHEMA_VERSION
    assert ThreeCompartmentRecovery.model_fields["schema_version"].default == SCHEMA_VERSION


def test_committed_schemas_match_the_models(repo_root) -> None:
    for name, model in CONTRACTS.items():
        committed = json.loads((repo_root / "schemas" / name).read_text(encoding="utf-8"))
        assert committed == json_schema(model)


def test_cross_project_fixtures_match_the_evidence_schema(repo_root) -> None:
    schema = json.loads((repo_root / "schemas" / "evidence_exchange.schema.json").read_text(encoding="utf-8"))
    for path in sorted((repo_root / "examples" / "cross_project_demo").glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        jsonschema.validate(instance=payload, schema=schema, format_checker=jsonschema.FormatChecker())
