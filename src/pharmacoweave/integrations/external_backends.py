"""Out-of-process engine boundary.

The only built-in executable backend is ``pharmacoweave-analytical``, a
MIT-licensed subprocess around the analytical exposure simulator. nlmixr2 and
PK-Sim stay in the catalog so their copyleft licenses are recorded, and neither
is imported or launched. A separate slot can start one caller-supplied
executable after that file, a version string, and license text are all present.
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
    CallerEngineRun,
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


_MAX_RESPONSE_BYTES = 2_000_000
_CALLER_LIMITATIONS = (
    "Research use only. The child response is stored as returned and is not a PharmacoWeave estimate.",
    "The license text was supplied by the caller and was not reviewed by this package.",
    "This slot does not import or vendor nlmixr2, PK-Sim, or NONMEM.",
    "pharmacoweave-analytical remains the only built-in executable backend.",
    "The run is not a dosing recommendation.",
)


def run_caller_engine(
    executable: str | Path,
    request_text: str,
    *,
    engine_version: str,
    license_text: str,
    timeout_s: float = 60.0,
) -> CallerEngineRun:
    """Start one caller-supplied executable when its file, version, and license text exist.

    The child receives the request path and the response path as its two arguments.
    A missing file, a blank version, or a blank license raises before any subprocess.
    """

    version = engine_version.strip()
    if not version or not license_text.strip():
        raise AnalysisError(
            "A caller engine starts only when the executable file, version, and license text are all present."
        )
    selected = Path(executable)
    if not selected.is_file() or selected.stat().st_size <= 0:
        raise AnalysisError(
            "A caller engine starts only when the executable file, version, and license text are all present."
        )
    file_bytes = selected.read_bytes()
    file_hash = hashlib.sha256(file_bytes).hexdigest()
    request_bytes = request_text.encode("utf-8")
    request_hash = hashlib.sha256(request_bytes).hexdigest()
    try:
        with tempfile.TemporaryDirectory(prefix="pharmacoweave-caller-engine-") as folder:
            root = Path(folder)
            request_path = root / "caller_request.txt"
            response_path = root / "caller_response.txt"
            request_path.write_bytes(request_bytes)
            completed = subprocess.run(
                [str(selected), str(request_path), str(response_path)],
                check=False,
                capture_output=True,
                timeout=timeout_s,
            )
            if completed.returncode != 0 or not response_path.is_file():
                detail = (completed.stderr or completed.stdout or b"no caller engine response").decode(
                    "utf-8", errors="replace"
                )
                raise AnalysisError(f"Caller engine exited {completed.returncode}: {detail.strip()[:2000]}")
            raw_response = response_path.read_bytes()
    except subprocess.TimeoutExpired as exc:
        raise AnalysisError(f"Caller engine exceeded {timeout_s} seconds.") from exc
    if len(raw_response) > _MAX_RESPONSE_BYTES:
        raise AnalysisError("Caller engine response exceeds the stored text limit.")
    try:
        response_text = raw_response.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AnalysisError("Caller engine response is not UTF-8 text.") from exc
    return CallerEngineRun(
        executable=str(selected.resolve()),
        engine_version=version,
        license_text=license_text,
        executable_sha256=file_hash,
        request_sha256=request_hash,
        response_sha256=hashlib.sha256(raw_response).hexdigest(),
        response_text=response_text,
        limitations=list(_CALLER_LIMITATIONS),
    )
