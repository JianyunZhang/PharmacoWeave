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
        "Single-dose IV bolus rich-sampling NCA with linear-up/log-down and an adjusted-R² terminal phase.",
    ),
    ModuleStatus(
        "pk_models",
        "P1",
        "partial",
        "IV bolus concentrations, log-linear OLS, and weighted NLS for CL and V. Other structures are not estimated.",
    ),
    ModuleStatus(
        "diagnostics",
        "P0",
        "partial",
        "Pass requires limitations. IV bolus fits write observed-vs-predicted and residual-vs-time SVG figures.",
    ),
    ModuleStatus(
        "reporting",
        "P0",
        "partial",
        "HTML validation report. NCA rows and IV bolus residual figures are included only after that analysis runs.",
    ),
    ModuleStatus(
        "integrations",
        "P0",
        "partial",
        "Evidence guards and a PharmacologyEvidence round trip with derived_from. No live Weave dependency.",
    ),
    ModuleStatus(
        "pd_models",
        "P1",
        "partial",
        "Direct-effect Emax with fixed baseline and hill 1. Sigmoid hill and indirect response are not estimated.",
    ),
    ModuleStatus(
        "pkpd",
        "P1",
        "partial",
        "Links an observed effect to an explicit concentration. AUC, Cmax, lag, and indirect response are not fitted.",
    ),
    ModuleStatus(
        "dmpk",
        "P1",
        "partial",
        "In vitro permeability, microsomal clearance, and hepatocyte clearance records. No in vivo scaling.",
    ),
    ModuleStatus(
        "toxicology",
        "P1",
        "partial",
        "Research-use tox endpoint records linked to a dose or exposure id. No clinical dose is produced.",
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
