"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.2.10 validates study data, calculates equal-dose interval NCA, and
calculates urine excretion amount and renal clearance from recorded void
volumes. A caller-supplied executable starts only when its file, version, and
license text are all present. A one-compartment constant-rate IV infusion can
be fit for CL and V. Finished individual estimates can be summarized on the
log scale, with or without inverse-variance weights from the caller standard error.
A caller physiology table can be checked, one linear turnover
equation has a closed form, and liver plus one tissue has an algebraic steady
state. In-process NLME, dynamic whole-body PBPK, and QSP networks stay refused.
A caller may supply a log-parameter
correlation for the prediction band. Steady-state clearance is Dose/AUC_tau
only when that steady state is declared. It does not give clinical dosing advice.
The only executable external backend is the MIT analytical reference engine.
nlmixr2 and PK-Sim are catalogued as copyleft and are not imported or launched.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
