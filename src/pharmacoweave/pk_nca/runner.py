"""NCA entry point.

The config and the result schema are executable checks. The numerical
calculator is intentionally absent until a reference terminal-phase rule is
implemented and tested.
"""

from pathlib import Path

import yaml
from pydantic import ValidationError

from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import IngestError, PlannedCapabilityError, ValidationFailed
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


def run_nca(study: PharmacologyStudy, config: str | Path) -> None:
    """Validate the study and config, then stop before calculating NCA."""

    loaded = load_nca_config(config)
    report = validate_study(study)
    if not report.ok:
        raise ValidationFailed(report)
    raise PlannedCapabilityError(
        module="pk_nca",
        phase="P0",
        detail=(
            f"Config accepted for {loaded.dosing} {loaded.design_route} {loaded.sampling} "
            f"with AUC method {loaded.auc_method}. The calculator is not implemented. "
            "Sparse sampling, urine, steady state, infusion, extravascular absorption, "
            "and multiple dosing are outside this config."
        ),
    )
