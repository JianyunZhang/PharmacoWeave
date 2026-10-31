"""PharmacoWeave, a research-use quantitative pharmacology toolkit.

Version 0.1.0 validates study data, units, explicit molecular-weight conversions,
auditable BLQ/replicate policies, and cross-project evidence boundaries. It
calculates single-dose rich-serial NCA for IV bolus, constant-rate IV infusion,
and extravascular (oral, subcutaneous, intramuscular) designs; estimates
1-compartment IV bolus, 1-compartment extravascular absorption, and
2-compartment IV bolus models with asymptotic standard errors and 95% confidence
intervals; fits direct and sigmoid Emax exposure-response curves on the original
effect scale; records plasma/microsomal unbound fractions and blood-to-plasma
ratios; runs exploratory well-stirred IVIVE only under explicit physiological
scaling parameters (output stays ``predicted``); audits prior predictions
against posterior in vivo estimates with fold-error metrics; bridges study
observations to typed effect and toxicology endpoints; aligns multi-scale
electrophysiology and omics endpoints with strict clock sync verification;
supports declarative raw table mapping (ingest_mapping.yaml) preserving source
file hashes; and verifies and compares ResultBundles offline. It does not run
sparse/multiple-dose designs, ship hidden species physiology tables, or give
clinical dosing advice.
"""

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__

__all__ = ["MODULE_STATUS", "PharmacologyStudy", "__version__"]
