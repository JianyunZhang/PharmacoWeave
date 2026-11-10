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
        "Study entities, closed units, explicit MW conversion, replicate collapsing, validation, "
        "observation bridge, bundles, verify-bundle, and diff.",
    ),
    ModuleStatus(
        "ingestion",
        "P0",
        "available",
        "Read and write study directories (CSV, TSV, Parquet) and declarative ingest_mapping.yaml "
        "preserving raw file SHA-256.",
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
        "Single-dose rich-serial NCA for IV bolus, IV infusion, and extravascular routes, "
        "equal-dose equal-interval and declared steady-state interval NCA, "
        "and urine excretion amount plus renal clearance kept separate from plasma AUC.",
    ),
    ModuleStatus(
        "pk_models",
        "P1",
        "partial",
        "1c IV bolus (OLS/WNLS), 1c IV infusion (CL, V), 1c extravascular (ka, CL/F, V/F), "
        "and 2c IV bolus (CL, V1, Q, V2) with SE and 95% CI.",
    ),
    ModuleStatus(
        "diagnostics",
        "P0",
        "partial",
        "Pass requires limitations. SVG diagnostics, plus prior-vs-estimated fold-error audit (AFE, AAFE, 2/3-fold).",
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
        "available",
        "Evidence guards, PharmacologyEvidence export, prediction audit, and multi-scale Ephys/Omics "
        "time alignment with clock sync checks. A caller-supplied executable starts only with a file, "
        "version, and license text; nlmixr2 and PK-Sim stay unlaunched.",
    ),
    ModuleStatus(
        "pd_models",
        "P1",
        "partial",
        "Original-scale nonlinear LS for direct and sigmoid Emax; study observation bridging to "
        "EffectObservation. Indirect response not fitted.",
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
        "In vitro permeability, microsomal/hepatocyte CLint, fu_plasma, fu_mic, and B:P records. "
        "Explicit-assumption well-stirred IVIVE (predicted grade only). No dose advice.",
    ),
    ModuleStatus(
        "toxicology",
        "P1",
        "partial",
        "Research-use tox endpoint records, schema, and study observation bridge. No clinical dose is produced.",
    ),
    ModuleStatus(
        "population_pk",
        "P2",
        "partial",
        "Two-stage geometric means of finished individual estimates, plus a reader for caller-supplied "
        "population tables. In-process NLME is refused. nlmixr2 is not imported or launched.",
    ),
    ModuleStatus(
        "pbpk",
        "P2",
        "partial",
        "Caller-supplied physiology checks and algebraic steady state for liver plus one tissue. "
        "Dynamic whole-body PBPK is refused. PK-Sim is not imported or launched.",
    ),
    ModuleStatus(
        "qsp",
        "P3",
        "partial",
        "Closed form for one linear turnover equation, evidence grade predicted. QSP networks are not simulated.",
    ),
    ModuleStatus(
        "simulation",
        "P4",
        "partial",
        "Single-dose, repeated-dose, and steady-state analytical exposure intervals. "
        "Optional caller-supplied log-parameter correlation. Not a dosing regimen or a population model.",
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
