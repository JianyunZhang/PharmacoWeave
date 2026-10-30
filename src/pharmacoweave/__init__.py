"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.0.4 validates study data, units, and cross-project evidence
boundaries. It calculates single-dose IV bolus NCA and estimates CL and V
for a one-compartment IV bolus by log-linear OLS or weighted nonlinear least
squares. It also records direct-effect Emax, in vitro ADME, and research-use
toxicology endpoints, and exports PharmacologyEvidence with derived_from.
It does not run other study designs, scale in vitro clearance in vivo, or
give dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
