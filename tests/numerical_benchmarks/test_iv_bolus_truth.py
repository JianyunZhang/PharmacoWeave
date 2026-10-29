"""Truth values the first NCA calculator has to recover.

NCA is not calculated here. The noise-free IV bolus fixture uses
CL = 0.20 L/h/kg and V = 1.0 L/kg, so a future rich-sampling NCA on the
quantified points, with the 24 h BLQ left censored, is expected to approach
AUC_inf = Dose / CL and t1/2 = ln(2) * V / CL.
"""

import math

import pytest
from examples.iv_bolus_case import LLOQ_UG_PER_ML, expected_concentration_ug_ml

from pharmacoweave.pk_models import auc_inf_iv_bolus, half_life


def test_subject_one_exposure_truth() -> None:
    body_kg = 0.25
    dose_mg = body_kg * 1.0
    clearance_l_h = body_kg * 0.20
    volume_l = body_kg * 1.0
    assert auc_inf_iv_bolus(
        dose=dose_mg,
        dose_unit="mg",
        clearance=clearance_l_h,
        clearance_unit="L/h",
        auc_unit="ug*h/mL",
    ) == pytest.approx(5.0)
    assert half_life(
        volume=volume_l,
        volume_unit="L",
        clearance=clearance_l_h,
        clearance_unit="L/h",
        time_unit="h",
    ) == pytest.approx(math.log(2.0) / 0.2)
    assert expected_concentration_ug_ml(body_kg, 720) > LLOQ_UG_PER_ML
    assert expected_concentration_ug_ml(body_kg, 1440) < LLOQ_UG_PER_ML
