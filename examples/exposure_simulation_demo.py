"""Exploratory single-dose exposure simulation.

The numbers are a declared one-compartment IV bolus, not an animal experiment
and not a dosing regimen.
"""

from pathlib import Path

from pharmacoweave.simulation.exposure import ExposureSimulationRequest, simulate_exposure


def main() -> None:
    request_path = Path(__file__).resolve().parent / "data" / "exposure_simulation_iv_bolus.json"
    request = ExposureSimulationRequest.model_validate_json(request_path.read_text(encoding="utf-8"))
    result = simulate_exposure(request)
    print(f"simulation_id={result.simulation_id}")
    print(f"Cmax={result.cmax} {result.concentration_unit} at {result.tmax_h} h")
    print(f"AUC_inf={result.auc_inf} {result.auc_unit}")
    print(f"half_life={result.half_life_h} h")
    print(f"time_above={result.time_above_threshold_h} h")
    print(f"evidence_type={result.evidence_type}")
    print(f"constitutes_dosing_recommendation={result.constitutes_dosing_recommendation}")


if __name__ == "__main__":
    main()
