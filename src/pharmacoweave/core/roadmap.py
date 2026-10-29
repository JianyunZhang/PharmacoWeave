"""Execution status of each package.

Reserved modules are importable so the repository layout matches the project
plan. Importing one does not run an analysis.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ModuleStatus:
    """Whether a package can be executed in this version."""

    name: str
    phase: str
    execution: str
    summary: str


MODULE_STATUS: tuple[ModuleStatus, ...] = (
    ModuleStatus(
        "core",
        "P0",
        "available",
        "Study entities, closed unit vocabulary, relational validation, and result bundles.",
    ),
    ModuleStatus(
        "ingestion",
        "P0",
        "available",
        "Read and write study directories as CSV, TSV, or Parquet plus study_metadata.json.",
    ),
    ModuleStatus(
        "bioanalysis",
        "P0",
        "available",
        "Auditable BLQ policies. The original observation status and number are preserved.",
    ),
    ModuleStatus(
        "pk_nca",
        "P1",
        "available",
        "Single-dose IV bolus rich-sampling NCA with linear-up/log-down. Other designs are rejected.",
    ),
    ModuleStatus(
        "pk_models",
        "P1",
        "available",
        "One-compartment IV bolus concentrations and CL/V estimation. Other structures are not estimated.",
    ),
    ModuleStatus(
        "diagnostics",
        "P1",
        "available",
        "IV bolus concentration-time, goodness-of-fit, and residual figures. Pass still requires limitations.",
    ),
    ModuleStatus(
        "reporting",
        "P1",
        "available",
        "Validation inventory, plus a research-use analysis report after NCA and the IV bolus fit run.",
    ),
    ModuleStatus(
        "integrations",
        "P0",
        "partial",
        "File-level evidence guards for OmicsWeave, TheraWeave, and EphysWeave. No live dependency.",
    ),
    ModuleStatus(
        "pd_models",
        "P1",
        "reserved",
        "Direct-effect, sigmoid Emax, and indirect-response models.",
    ),
    ModuleStatus(
        "pkpd",
        "P1",
        "reserved",
        "Joint concentration-effect estimation, lag, and indirect-response linkage.",
    ),
    ModuleStatus(
        "dmpk",
        "P1",
        "reserved",
        "In vitro stability, binding, clearance, and tissue-distribution results.",
    ),
    ModuleStatus(
        "toxicology",
        "P1",
        "reserved",
        "Toxicology endpoints linked to dose and exposure, including TK.",
    ),
    ModuleStatus(
        "population_pk",
        "P2",
        "reserved",
        "Population PK through an external NLME engine, not an in-house NONMEM replacement.",
    ),
    ModuleStatus(
        "pbpk",
        "P2",
        "reserved",
        "PBPK through an external engine adapter.",
    ),
    ModuleStatus(
        "qsp",
        "P3",
        "reserved",
        "Quantitative systems pharmacology models.",
    ),
    ModuleStatus(
        "simulation",
        "P2",
        "reserved",
        "Study-design simulation under an explicit context of use.",
    ),
)

_BY_NAME = {item.name: item for item in MODULE_STATUS}


def status(name: str) -> ModuleStatus:
    """Return the execution status recorded for ``name``."""

    try:
        return _BY_NAME[name]
    except KeyError as exc:
        known = ", ".join(sorted(_BY_NAME))
        raise KeyError(f"Unknown module {name!r}. Known modules: {known}") from exc
