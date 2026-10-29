"""NCA on the synthetic IV bolus study, against a frozen trapezoidal reference.

Tolerances are also written in ``tests/numerical_benchmarks/README.md`` and
``docs/validation.md``. The closed-form Dose/CL value is a separate check. It
is not what the calculator returns.
"""

import importlib.util
import math
from pathlib import Path

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.pk_nca import run_nca

_SPEC = importlib.util.spec_from_file_location(
    "reference_linlog",
    Path(__file__).with_name("reference_linlog.py"),
)
assert _SPEC is not None and _SPEC.loader is not None
_REFERENCE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_REFERENCE)
reference_nca = _REFERENCE.reference_nca

# Frozen output of reference_linlog.reference_nca on the 6-decimal quantified
# curve. Regenerating these from the library would not be an external check.
FROZEN = {
    "auc_last": 4.546412284265698,
    "auc_inf": 5.000002593636931,
    "fraction": 0.09071801481632784,
    "cmax": 0.983471,
    "tmax_h": 5 / 60,
    "lambda_z": 0.19999986359001665,
    "half_life_h": 3.465738266606223,
    "n_points": 7,
    "adjusted_r2": 0.99999999999912,
}
_TIMES_H = [minute / 60 for minute in (5, 15, 30, 60, 120, 240, 480, 720)]
_CONCS = [0.983471, 0.951229, 0.904837, 0.818731, 0.670320, 0.449329, 0.201897, 0.090718]
# Absolute acceptance versus the noise-free analytical targets. The gap is the
# 6-decimal rounding of the fixture, not a substitution of Dose/CL for AUC.
_ANALYTICAL_AUC_ABS = 1e-4
_ANALYTICAL_LAMBDA_ABS = 1e-5
_ANALYTICAL_HALF_LIFE_ABS = 1e-4
_REFERENCE_ABS = 1e-9


def test_frozen_reference_matches_the_independent_transcription() -> None:
    produced = reference_nca(_TIMES_H, _CONCS)
    for key, expected in FROZEN.items():
        assert produced[key] == pytest.approx(expected, abs=_REFERENCE_ABS)


def test_library_matches_frozen_trapezoids_and_analytical_targets(repo_root) -> None:
    study = build_study()
    before = [(item.observation_id, item.status, item.numeric_value) for item in study.observations]
    records = run_nca(study, repo_root / "configs" / "nca_single_dose_iv.yaml")
    after = [(item.observation_id, item.status, item.numeric_value) for item in study.observations]
    assert after == before
    assert len(records.nca_parameters) == 2
    for record, dose_mg, clearance_l_h, volume_l in (
        (records.nca_parameters[0], 0.25, 0.05, 0.25),
        (records.nca_parameters[1], 0.3, 0.06, 0.30),
    ):
        assert record.terminal_phase_qualified is True
        assert record.auc_last == pytest.approx(FROZEN["auc_last"], abs=_REFERENCE_ABS)
        assert record.auc_inf == pytest.approx(FROZEN["auc_inf"], abs=_REFERENCE_ABS)
        assert record.auc_inf != 5.0
        assert record.auc_inf == pytest.approx(5.0, abs=_ANALYTICAL_AUC_ABS)
        assert record.auc_inf_extrap_fraction == pytest.approx(FROZEN["fraction"], abs=_REFERENCE_ABS)
        assert record.cmax == pytest.approx(FROZEN["cmax"], abs=_REFERENCE_ABS)
        assert record.tmax == pytest.approx(FROZEN["tmax_h"], abs=_REFERENCE_ABS)
        assert record.lambda_z == pytest.approx(FROZEN["lambda_z"], abs=_REFERENCE_ABS)
        assert record.lambda_z == pytest.approx(0.2, abs=_ANALYTICAL_LAMBDA_ABS)
        assert record.half_life == pytest.approx(FROZEN["half_life_h"], abs=_REFERENCE_ABS)
        assert record.half_life == pytest.approx(math.log(2.0) / 0.2, abs=_ANALYTICAL_HALF_LIFE_ABS)
        assert record.terminal_phase_n_points == FROZEN["n_points"]
        assert record.clearance == pytest.approx(clearance_l_h, rel=1e-4)
        assert record.clearance == pytest.approx(dose_mg / FROZEN["auc_inf"], rel=1e-9)
        assert record.volume == pytest.approx(volume_l, rel=1e-4)
        assert record.clearance_semantic is not None and record.clearance_semantic.value == "CL"
        assert record.volume_semantic is not None and record.volume_semantic.value == "Vz"
        assert any("not set to zero" in note for note in record.warnings)
        zero_imputed = FROZEN["auc_last"] + 0.5 * _CONCS[-1] * (24 - 12)
        assert abs((record.auc_last or 0.0) - zero_imputed) > 0.4
    blq = [item for item in study.observations if item.sample_id.endswith("T1440")]
    assert len(blq) == 2
    assert all(item.status is ObservationStatus.BLQ and item.numeric_value is None for item in blq)
