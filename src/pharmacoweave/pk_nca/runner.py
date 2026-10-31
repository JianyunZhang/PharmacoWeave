"""NCA entry point for supported single-dose rich-serial designs."""

from pathlib import Path

import yaml
from pydantic import ValidationError

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.core.profiles import single_dose_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import IngestError, ValidationFailed
from pharmacoweave.pk_nca.calculate import nca_profile
from pharmacoweave.pk_nca.config import NcaConfig


def load_nca_config(path: str | Path) -> NcaConfig:
    """Load an NCA config and reject designs outside the locked scope."""

    selected = Path(path)
    if not selected.is_file():
        raise IngestError(f"NCA config {selected} does not exist.")
    try:
        payload = yaml.safe_load(selected.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise IngestError(f"{selected} is not valid YAML: {exc}") from exc
    if not isinstance(payload, dict):
        raise IngestError(f"{selected} must contain a mapping.")
    try:
        return NcaConfig.model_validate(payload)
    except ValidationError as exc:
        raise IngestError(f"{selected} is outside the implemented NCA design: {exc}") from exc


def run_nca(study: PharmacologyStudy, config: str | Path | NcaConfig) -> list[NcaParameterSet]:
    """Run single-dose rich-serial NCA for the configured route.

    Sparse sampling, urine, steady state, and multiple dosing are rejected by
    the config or by the profile builder.
    """

    loaded = config if isinstance(config, NcaConfig) else load_nca_config(config)
    report = validate_study(study)
    if not report.ok:
        raise ValidationFailed(report)
    expected_route = Route(loaded.design_route)
    profiles = single_dose_profiles(study, allowed_routes=(expected_route,))
    return [nca_profile(profile, loaded) for profile in profiles]
