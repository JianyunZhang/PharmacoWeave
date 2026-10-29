"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.0.1 validates study data, units, and cross-project evidence
boundaries. It calculates one-compartment IV bolus concentrations from known
parameters. It does not estimate NCA parameters or fit PK models.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
