"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.1.0 validates study data and calculates single-dose IV bolus NCA
plus a one-compartment fit for CL and V. It does not produce dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
