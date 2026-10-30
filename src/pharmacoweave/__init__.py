"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.0.2 validates study data, units, and cross-project evidence
boundaries. It calculates single-dose IV bolus NCA and a log-linear
one-compartment fit for CL and V. It does not run other study designs.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
