"""NCA entry point for the locked single-dose IV bolus design.

Other designs are rejected. The calculator does not run for sparse sampling,
urine, steady state, infusion, extravascular absorption, or multiple dosing.
"""

from pathlib import Path

import yaml
from pydantic import ValidationError

from pharmacoweave.core.modeling import AnalysisRecords
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import IngestError, ValidationFailed
from pharmacoweave.pk_nca.calculator import calculate_profile
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_profiles import load_iv_bolus_profiles


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


def run_nca(study: PharmacologyStudy, config: NcaConfig | str | Path) -> AnalysisRecords:
    """Calculate NCA for a validated single-dose IV bolus study.

    The returned records contain one parameter set per concentration profile.
    Terminal parameters and AUC∞ are absent when their reporting rules fail.
    BLQ rows on the study are not modified.
    """

    loaded = config if isinstance(config, NcaConfig) else load_nca_config(config)
    report = validate_study(study)
    if not report.ok:
        raise ValidationFailed(report)
    profiles = load_iv_bolus_profiles(study)
    return AnalysisRecords(nca_parameters=[calculate_profile(profile, loaded) for profile in profiles])
