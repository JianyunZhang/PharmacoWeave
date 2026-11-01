"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.2.0 validates study data, units, explicit molecular-weight conversions,
auditable BLQ/replicate policies, and cross-project evidence boundaries. It
calculates single-dose rich-serial NCA for IV bolus, constant-rate IV infusion,
and extravascular designs; estimates 1-compartment IV, 1-compartment
extravascular, and 2-compartment IV models; fits direct and sigmoid Emax curves;
runs exploratory well-stirred IVIVE and prediction-audit; and simulates
single-dose analytical exposure intervals under an explicit context of use.
The only executable external backend is the MIT analytical reference engine.
nlmixr2 and PK-Sim are catalogued as copyleft and are not imported or launched.
It does not give clinical dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
