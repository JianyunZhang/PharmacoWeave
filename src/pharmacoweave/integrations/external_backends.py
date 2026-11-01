"""Out-of-process engine boundary.

The only executable backend is ``pharmacoweave-analytical``, a MIT-licensed
subprocess around the analytical exposure simulator. nlmixr2 and PK-Sim stay
in the catalog so their copyleft licenses are recorded, and neither is imported
or launched.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from pharmacoweave.errors import AnalysisError, PlannedCapabilityError
from pharmacoweave.simulation.engine_protocol import (
    ANALYTICAL_ENGINE,
    ENGINE_CATALOG,
    EngineCatalogEntry,
    EngineRequest,
    EngineResponse,
    catalog_entry,
)
from pharmacoweave.simulation.exposure import ExposureSimulationRequest


def supported_backends() -> list[str]:
    """Return engine names this package can execute."""

    return [entry.name for entry in ENGINE_CATALOG if entry.executable]


def describe_engine(name: str) -> EngineCatalogEntry:
    """Return the license and execution flag for a known engine without starting it."""

    return catalog_entry(name)


def run_backend(
    name: str,
    request: ExposureSimulationRequest | None = None,
    *,
    timeout_s: float = 60.0,
) -> EngineResponse:
    """Run an executable backend out of process, or refuse a copyleft engine."""

    entry = catalog_entry(name)
    if not entry.executable:
        raise PlannedCapabilityError(
            module="integrations.external_backends",
            phase="P4",
            detail=(
                f"{name} is registered as {entry.license_spdx} (copyleft={entry.copyleft}) and is not executable. "
                f"{entry.limitation} Supported backends: {', '.join(supported_backends())}."
            ),
        )
    if name != ANALYTICAL_ENGINE:
        raise PlannedCapabilityError(
            module="integrations.external_backends",
            phase="P4",
            detail=f"{name} has no subprocess adapter in this version.",
        )
    if request is None:
        raise AnalysisError(f"{name} requires an ExposureSimulationRequest.")
    payload = EngineRequest(simulation=request).model_dump_json()
    src = str(Path(__file__).resolve().parents[2])
    env = os.environ.copy()
    previous = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = src if not previous else src + os.pathsep + previous
    with tempfile.TemporaryDirectory(prefix="pharmacoweave-engine-") as folder:
        root = Path(folder)
        request_path = root / "engine_request.json"
        response_path = root / "engine_response.json"
        request_path.write_text(payload, encoding="utf-8")
        completed = subprocess.run(
            [sys.executable, "-m", "pharmacoweave.simulation.reference_engine", str(request_path), str(response_path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            env=env,
        )
        if completed.returncode != 0 or not response_path.is_file():
            detail = (completed.stderr or completed.stdout or "no engine response").strip()
            raise AnalysisError(f"{name} exited {completed.returncode}: {detail[:2000]}")
        response = EngineResponse.model_validate_json(response_path.read_text(encoding="utf-8"))
    expected_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    if (
        response.engine_name != name
        or response.copyleft
        or response.license_spdx != "MIT"
        or response.constitutes_dosing_recommendation
        or response.request_sha256 != expected_hash
    ):
        raise AnalysisError("Engine response failed the PharmacoWeave exchange contract.")
    return response
