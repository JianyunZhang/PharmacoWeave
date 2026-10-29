"""Checked-in JSON Schemas stay aligned with the Pydantic models."""

import json

import jsonschema

from pharmacoweave.core.provenance import ResultBundleManifest
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.schema_export import CONTRACTS, json_schema
from pharmacoweave.version import SCHEMA_VERSION


def test_schema_version_defaults_match_the_package_constant() -> None:
    assert ResultBundleManifest.model_fields["schema_version"].default == SCHEMA_VERSION
    assert NcaConfig.model_fields["schema_version"].default == SCHEMA_VERSION


def test_committed_schemas_match_the_models(repo_root) -> None:
    for name, model in CONTRACTS.items():
        committed = json.loads((repo_root / "schemas" / name).read_text(encoding="utf-8"))
        assert committed == json_schema(model)


def test_cross_project_fixtures_match_the_evidence_schema(repo_root) -> None:
    schema = json.loads((repo_root / "schemas" / "evidence_exchange.schema.json").read_text(encoding="utf-8"))
    for path in sorted((repo_root / "examples" / "cross_project_demo").glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        jsonschema.validate(instance=payload, schema=schema, format_checker=jsonschema.FormatChecker())
