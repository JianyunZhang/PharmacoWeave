"""Single-dose rich-serial noncompartmental analysis."""

from pharmacoweave.core.modeling import NcaParameterSet, auc_inf_is_reportable
from pharmacoweave.core.roadmap import status
from pharmacoweave.pk_nca.calculate import (
    nca_extravascular,
    nca_iv_bolus,
    nca_iv_infusion,
    nca_profile,
)
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_nca.runner import load_nca_config, run_nca

STATUS = status("pk_nca")

__all__ = [
    "NcaConfig",
    "NcaParameterSet",
    "STATUS",
    "auc_inf_is_reportable",
    "load_nca_config",
    "nca_extravascular",
    "nca_iv_bolus",
    "nca_iv_infusion",
    "nca_profile",
    "run_nca",
]
