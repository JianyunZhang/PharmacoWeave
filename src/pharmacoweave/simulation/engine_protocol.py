"""File-exchange records for an out-of-process analytical reference engine.

nlmixr2 and PK-Sim are described here so a result can cite their license.
They are not imported, and this module does not start them.
"""

import hashlib
from typing import Literal

from pharmacoweave.core.entities import Record
from pharmacoweave.errors import AnalysisError
from pharmacoweave.simulation.exposure import ExposureSimulationRequest, ExposureSimulationResult
from pharmacoweave.version import __version__

ANALYTICAL_ENGINE: Literal["pharmacoweave-analytical"] = "pharmacoweave-analytical"


class EngineCatalogEntry(Record):
    """A known engine. ``executable`` is true only for the in-repo reference engine."""

    name: str
    license_spdx: str
    copyleft: bool
    executable: bool
    version_note: str
    limitation: str


class EngineRequest(Record):
    """JSON file handed to the reference-engine subprocess."""

    schema_name: Literal["pharmacoweave-engine-request"] = "pharmacoweave-engine-request"
    schema_version: Literal["0.2.1"] = "0.2.1"
    engine: Literal["pharmacoweave-analytical"] = ANALYTICAL_ENGINE
    simulation: ExposureSimulationRequest


class EngineResponse(Record):
    """JSON file returned by the reference-engine subprocess."""

    schema_name: Literal["pharmacoweave-engine-exchange"] = "pharmacoweave-engine-exchange"
    schema_version: Literal["0.2.1"] = "0.2.1"
    engine_name: Literal["pharmacoweave-analytical"] = ANALYTICAL_ENGINE
    engine_version: str
    license_spdx: Literal["MIT"] = "MIT"
    copyleft: Literal[False] = False
    executed: Literal[True] = True
    context_of_use: str
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    request_sha256: str
    result: ExposureSimulationResult


ENGINE_CATALOG: tuple[EngineCatalogEntry, ...] = (
    EngineCatalogEntry(
        name=ANALYTICAL_ENGINE,
        license_spdx="MIT",
        copyleft=False,
        executable=True,
        version_note=__version__,
        limitation=(
            "Subprocess file exchange around the analytical single-dose simulator already tested in this package. "
            "Not a population or PBPK engine."
        ),
    ),
    EngineCatalogEntry(
        name="nlmixr2",
        license_spdx="GPL-2.0-or-later",
        copyleft=True,
        executable=False,
        version_note="not-bundled",
        limitation=(
            "Copyleft NLME engine. PharmacoWeave does not import it and does not start a subprocess for it "
            "in this version."
        ),
    ),
    EngineCatalogEntry(
        name="pk-sim",
        license_spdx="GPL-2.0-only",
        copyleft=True,
        executable=False,
        version_note="not-bundled",
        limitation=(
            "Copyleft PBPK engine from Open Systems Pharmacology. PharmacoWeave does not import it and does not "
            "start a subprocess for it in this version."
        ),
    ),
)


def catalog_entry(name: str) -> EngineCatalogEntry:
    """Return the catalog row for ``name``."""

    for entry in ENGINE_CATALOG:
        if entry.name == name:
            return entry
    known = ", ".join(entry.name for entry in ENGINE_CATALOG)
    raise AnalysisError(f"Unknown engine {name!r}. Known engines: {known}.")


def response_for(raw_request: str, result: ExposureSimulationResult) -> EngineResponse:
    """Bind a simulation result to the exact request bytes the child received."""

    return EngineResponse(
        engine_version=__version__,
        context_of_use=result.context_of_use,
        limitations=list(result.limitations)
        + [
            f"Executed out of process by {ANALYTICAL_ENGINE} {__version__} (MIT, not copyleft).",
            "nlmixr2 and PK-Sim were not imported and were not started.",
        ],
        request_sha256=hashlib.sha256(raw_request.encode("utf-8")).hexdigest(),
        result=result,
    )
