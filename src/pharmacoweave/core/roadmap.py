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
        "Study entities, closed unit vocabulary, replicate collapsing, validation, and result bundles.",
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
        "Auditable BLQ policies and derived observation tables without mutating raw observations.",
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
        "IV bolus log-linear OLS and weighted NLS (with asymptotic SE and 95% CI) for CL and V.",
    ),
    ModuleStatus(
        "diagnostics",
        "P0",
        "partial",
        "Pass requires limitations. Writes semi-log, GoF, residual, weighted-residual, and Emax SVG diagnostics.",
    ),
    ModuleStatus(
        "reporting",
        "P0",
        "partial",
        "HTML report and native SVGs for NCA, PK compartment fits, and direct/sigmoid Emax curves.",
    ),
    ModuleStatus(
        "integrations",
        "P0",
        "partial",
        "Evidence guards and PharmacologyEvidence export across NCA, PK fits, Emax, ADME, and toxicology.",
    ),
    ModuleStatus(
        "pd_models",
        "P1",
        "partial",
        "Original-scale nonlinear LS for direct and sigmoid Emax (with SE and 95% CI). Indirect response not fitted.",
    ),
    ModuleStatus(
        "pkpd",
        "P1",
        "partial",
        "Links observed direct effects to explicit concentrations for direct/sigmoid Emax. Lag and indirect refused.",
    ),
    ModuleStatus(
        "dmpk",
        "P1",
        "partial",
        "In vitro permeability, microsomal clearance, and hepatocyte clearance records and schema. No IVIVE scaling.",
    ),
    ModuleStatus(
        "toxicology",
        "P1",
        "partial",
        "Research-use tox endpoint records and schema linked to dose or exposure id. No clinical dose is produced.",
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
