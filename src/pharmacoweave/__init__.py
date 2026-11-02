"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.2.1 validates study data and simulates single-dose, repeated-dose, and
steady-state analytical exposure intervals under an explicit context of use.
The only executable external backend is the MIT analytical reference engine.
nlmixr2 and PK-Sim are catalogued as copyleft and are not imported or launched.
It does not give clinical dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
