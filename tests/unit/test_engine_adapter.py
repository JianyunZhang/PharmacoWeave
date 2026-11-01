"""The executable engine is the MIT analytical subprocess. Copyleft engines stay refused."""

import ast
import subprocess
from pathlib import Path

import pytest

from pharmacoweave.errors import AnalysisError, PlannedCapabilityError
from pharmacoweave.integrations.external_backends import describe_engine, run_backend, supported_backends
from pharmacoweave.pbpk import run_in_process as run_pbpk
from pharmacoweave.population_pk import run_in_process as run_population
from pharmacoweave.qsp import run_in_process as run_qsp
from pharmacoweave.simulation.exposure import ExposureSimulationRequest, simulate_exposure


def _request() -> ExposureSimulationRequest:
    return ExposureSimulationRequest.model_validate(
        {
            "simulation_id": "SIM-ENGINE",
            "context_of_use": "Research check of the out-of-process analytical reference engine.",
            "model_id": "one_compartment_iv_bolus",
            "dose": 0.25,
            "dose_unit": "mg",
            "clearance": 0.05,
            "clearance_unit": "L/h",
            "volume": 0.05,
            "volume_unit": "L",
            "times_h": [0, 1, 4],
        }
    )


def test_only_the_analytical_engine_is_executable() -> None:
    assert supported_backends() == ["pharmacoweave-analytical"]
    nlmixr = describe_engine("nlmixr2")
    pksim = describe_engine("pk-sim")
    analytical = describe_engine("pharmacoweave-analytical")
    assert nlmixr.copyleft and nlmixr.executable is False
    assert pksim.copyleft and pksim.license_spdx == "GPL-2.0-only"
    assert analytical.executable and analytical.license_spdx == "MIT"
    with pytest.raises(PlannedCapabilityError, match="not executable"):
        run_backend("nlmixr2", _request())
    with pytest.raises(PlannedCapabilityError):
        run_backend("pk-sim")
    with pytest.raises(AnalysisError, match="Unknown engine"):
        describe_engine("nonmem")


def test_copyleft_engine_does_not_start_a_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("subprocess was started for a copyleft engine")

    monkeypatch.setattr(subprocess, "run", fail_if_called)
    with pytest.raises(PlannedCapabilityError):
        run_backend("nlmixr2", _request())


def test_reference_engine_matches_the_in_process_curve() -> None:
    request = _request()
    direct = simulate_exposure(request)
    response = run_backend("pharmacoweave-analytical", request)
    assert response.copyleft is False
    assert response.license_spdx == "MIT"
    assert response.constitutes_dosing_recommendation is False
    assert response.result.cmax == pytest.approx(direct.cmax, rel=1e-12)
    assert response.result.auc_inf == pytest.approx(direct.auc_inf, rel=1e-12)
    assert response.result.curve[1].concentration == pytest.approx(direct.curve[1].concentration, rel=1e-12)
    assert "nlmixr2 and PK-Sim were not imported" in " ".join(response.limitations)


def test_reserved_model_families_stay_refused() -> None:
    with pytest.raises(PlannedCapabilityError):
        run_population()
    with pytest.raises(PlannedCapabilityError):
        run_pbpk()
    with pytest.raises(PlannedCapabilityError):
        run_qsp()


def test_engine_modules_do_not_import_copyleft_solvers() -> None:
    forbidden = {"nlmixr2", "nlmixr", "rpy2", "pksim", "nonmem"}
    root = Path(__file__).resolve().parents[2] / "src" / "pharmacoweave"
    for relative in (
        "simulation/exposure.py",
        "simulation/engine_protocol.py",
        "simulation/reference_engine.py",
        "integrations/external_backends.py",
    ):
        tree = ast.parse((root / relative).read_text(encoding="utf-8"))
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0].lower() for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0].lower())
        assert imported.isdisjoint(forbidden)
