"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.0.5 validates study data, units, auditable BLQ/replicate policies,
and cross-project evidence boundaries. It calculates single-dose IV bolus NCA,
estimates CL and V for a one-compartment IV bolus by log-linear OLS or weighted
NLS with asymptotic standard errors and 95% confidence intervals, fits direct
and sigmoid Emax exposure-response curves on the original effect scale, and
persists Phase 2 PD/ADME/toxicology tables, diagnostic SVGs, and screened
PharmacologyEvidence inside checksummed result bundles. It does not run other
study designs, scale in vitro clearance in vivo, or give dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
