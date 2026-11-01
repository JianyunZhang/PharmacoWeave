"""Subprocess entry point for the analytical reference engine.

Run as ``python -m pharmacoweave.simulation.reference_engine REQUEST RESPONSE``.
The process reads one JSON request and writes one JSON response. It does not
import nlmixr2 or PK-Sim.
"""

from __future__ import annotations

import sys
from pathlib import Path

from pharmacoweave.simulation.engine_protocol import ANALYTICAL_ENGINE, EngineRequest, response_for
from pharmacoweave.simulation.exposure import simulate_exposure


def main(argv: list[str] | None = None) -> int:
    """Read an engine request path and write an engine response path."""

    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        print(
            "usage: python -m pharmacoweave.simulation.reference_engine REQUEST RESPONSE",
            file=sys.stderr,
        )
        return 2
    request_path = Path(args[0])
    response_path = Path(args[1])
    raw = request_path.read_text(encoding="utf-8")
    request = EngineRequest.model_validate_json(raw)
    if request.engine != ANALYTICAL_ENGINE:
        print(f"refusing engine {request.engine!r}", file=sys.stderr)
        return 2
    result = simulate_exposure(request.simulation)
    response_path.write_text(response_for(raw, result).model_dump_json(indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
