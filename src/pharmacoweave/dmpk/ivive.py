"""Exploratory, explicit-assumption well-stirred hepatic IVIVE.

The calculation follows the well-stirred liver model (e.g. Yang et al.,
2007, Curr Drug Metab 8:676-684):

- ``CLint,u,vitro = CLint,app / fu_inc``
- ``CLint,u,H [L/h] = CLint,u,vitro * (MPPGL or HPGL) * liver_weight_g * 60 / 1e6``
- ``fu_b = fu_p / (B:P)``
- ``CLH,b = QH * fu_b * CLint,u,H / (QH + fu_b * CLint,u,H)``
- ``EH = CLH,b / QH`` and ``CLH,p = (B:P) * CLH,b``

All physiological inputs come from the caller. A scalar ``fu_inc`` or
``blood_to_plasma_ratio`` is accepted only as an explicit assumption and is
written into ``assumptions``. The result is a ``predicted`` record. It never
becomes an observed or estimated CL, and no starting dose is derived from it.
"""

from typing import NoReturn

from pharmacoweave.core.entities import Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.dmpk.physiology import (
    IVIVE_EXPLORATORY_LIMITATION,
    HepaticClearanceModel,
    IviveClearancePrediction,
    PhysiologicalScalingParams,
)
from pharmacoweave.dmpk.records import AdmeAssayKind, AdmeStatus, InVitroAdmeRecord
from pharmacoweave.errors import AnalysisError, ContractViolation

_MIN_PER_H = 60.0
_UL_PER_L = 1e6


class IviveRequest(Record):
    """Input payload for exploratory well-stirred IVIVE calculation."""

    clint_record: InVitroAdmeRecord
    physiology: PhysiologicalScalingParams
    fu_plasma: float | InVitroAdmeRecord
    fu_inc: float | InVitroAdmeRecord = 1.0
    blood_to_plasma_ratio: float | InVitroAdmeRecord = 1.0
    prediction_id: str | None = None


def predict_hepatic_clearance_well_stirred(
    clint_record: InVitroAdmeRecord,
    *,
    physiology: PhysiologicalScalingParams,
    fu_plasma: InVitroAdmeRecord | float,
    fu_inc: InVitroAdmeRecord | float,
    blood_to_plasma_ratio: InVitroAdmeRecord | float,
    prediction_id: str | None = None,
) -> IviveClearancePrediction:
    """Return an exploratory well-stirred hepatic clearance prediction.

    ``clint_record`` must be a measured microsomal or hepatocyte clearance.
    ``fu_plasma`` is a measured ``fu_plasma`` record or an explicit fraction in
    (0, 1]. ``fu_inc`` is a measured ``fu_mic`` record (microsomes only) or an
    explicit fraction. ``blood_to_plasma_ratio`` is a measured record or an
    explicit positive ratio. Nothing is defaulted.
    """

    if clint_record.assay_kind not in {AdmeAssayKind.MICROSOMAL_CLEARANCE, AdmeAssayKind.HEPATOCYTE_CLEARANCE}:
        raise AnalysisError(
            f"{clint_record.record_id}: IVIVE needs microsomal_clearance or hepatocyte_clearance, "
            f"not {clint_record.assay_kind.value}."
        )
    clint = clint_record.canonical_value()
    if clint is None:
        raise AnalysisError(
            f"{clint_record.record_id} is {clint_record.status.value}; only a measured CLint can be scaled."
        )
    clint_value, clint_unit = clint
    microsomal = clint_record.assay_kind is AdmeAssayKind.MICROSOMAL_CLEARANCE
    if microsomal:
        scaling = physiology.microsomal_protein_mg_per_g_liver
        scaling_unit = "mg_protein/g_liver"
        if scaling is None:
            raise AnalysisError(
                f"Physiology profile {physiology.profile_id} has no microsomal_protein_mg_per_g_liver (MPPGL). "
                "It is not assumed."
            )
    else:
        scaling = physiology.hepatocellularity_million_cells_per_g_liver
        scaling_unit = "million_cells/g_liver"
        if scaling is None:
            raise AnalysisError(
                f"Physiology profile {physiology.profile_id} has no hepatocellularity_million_cells_per_g_liver "
                "(HPGL). It is not assumed."
            )

    assumptions: list[str] = [
        "Well-stirred liver model; hepatic metabolism is the only clearance pathway scaled.",
        "Linear (non-saturated) kinetics at the in vitro substrate concentration.",
        (
            f"Physiology profile {physiology.profile_id} ({physiology.species}): body weight "
            f"{physiology.body_weight_kg:g} kg, liver {physiology.liver_weight_g:g} g, QH "
            f"{physiology.hepatic_blood_flow_l_h:g} L/h, {scaling_unit.replace('_', ' ')} {scaling:g}. "
            f"Source: {physiology.source}"
        ),
        f"In vitro system: {clint_record.system}.",
    ]
    derived_from = [
        f"adme:{clint_record.record_id}",
        f"compound:{clint_record.compound_id}",
        f"physiology:{physiology.profile_id}",
    ]

    fu_p = _fraction_input(
        fu_plasma,
        label="fu_plasma",
        expected_kind=AdmeAssayKind.FU_PLASMA,
        compound_id=clint_record.compound_id,
        assumptions=assumptions,
        derived_from=derived_from,
    )
    if isinstance(fu_inc, InVitroAdmeRecord) and not microsomal:
        raise AnalysisError(
            "A fu_mic record describes microsomal binding and is not used for a hepatocyte incubation. "
            "Pass fu_inc for hepatocytes as an explicit scalar assumption."
        )
    fu_i = _fraction_input(
        fu_inc,
        label="fu_inc",
        expected_kind=AdmeAssayKind.FU_MIC,
        compound_id=clint_record.compound_id,
        assumptions=assumptions,
        derived_from=derived_from,
    )
    ratio = _ratio_input(
        blood_to_plasma_ratio,
        compound_id=clint_record.compound_id,
        assumptions=assumptions,
        derived_from=derived_from,
    )

    qh = physiology.hepatic_blood_flow_l_h
    clint_u_vitro = clint_value / fu_i
    clint_u_h = clint_u_vitro * scaling * physiology.liver_weight_g * _MIN_PER_H / _UL_PER_L
    fu_b = fu_p / ratio
    cl_blood = qh * fu_b * clint_u_h / (qh + fu_b * clint_u_h)
    extraction = cl_blood / qh
    cl_plasma = ratio * cl_blood
    if not 0.0 < extraction < 1.0:
        raise AnalysisError(f"Hepatic extraction ratio {extraction} is outside (0, 1); check the inputs.")

    return IviveClearancePrediction(
        prediction_id=prediction_id or f"IVIVE-{clint_record.record_id}-{physiology.profile_id}",
        compound_id=clint_record.compound_id,
        species=physiology.species,
        physiology_profile_id=physiology.profile_id,
        source_record_id=clint_record.record_id,
        assay_kind="microsomal_clearance" if microsomal else "hepatocyte_clearance",
        clint_app_value=clint_value,
        clint_app_unit="uL/min/mg_protein" if clint_unit == "uL/min/mg_protein" else "uL/min/million_cells",
        fu_inc=fu_i,
        fu_plasma=fu_p,
        blood_to_plasma_ratio=ratio,
        fu_blood=fu_b,
        clint_u_vitro_value=clint_u_vitro,
        scaling_factor_value=scaling,
        scaling_factor_unit="mg_protein/g_liver" if microsomal else "million_cells/g_liver",
        liver_weight_g=physiology.liver_weight_g,
        body_weight_kg=physiology.body_weight_kg,
        hepatic_blood_flow_l_h=qh,
        clint_u_hepatic_l_h=clint_u_h,
        hepatic_extraction_ratio=extraction,
        hepatic_blood_clearance_l_h=cl_blood,
        hepatic_plasma_clearance_l_h=cl_plasma,
        hepatic_plasma_clearance_l_h_kg=cl_plasma / physiology.body_weight_kg,
        model=HepaticClearanceModel.WELL_STIRRED,
        assumptions=assumptions,
        limitations=[
            RESEARCH_USE_LIMITATION,
            IVIVE_EXPLORATORY_LIMITATION,
            "Hepatic clearance only. Renal, biliary, and extrahepatic clearance are not included.",
        ],
        derived_from=derived_from,
    )


def human_starting_dose_from_ivive(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse to turn an IVIVE prediction into a starting dose or FIH recommendation."""

    raise ContractViolation(
        "PharmacoWeave does not derive a human starting dose, MRSD, or FIH dose from an IVIVE prediction. "
        "IVIVE output stays an exploratory predicted record."
    )


def _fraction_input(
    value: InVitroAdmeRecord | float,
    *,
    label: str,
    expected_kind: AdmeAssayKind,
    compound_id: str,
    assumptions: list[str],
    derived_from: list[str],
) -> float:
    if isinstance(value, InVitroAdmeRecord):
        fraction = _measured_record_value(value, expected_kind=expected_kind, compound_id=compound_id)
        derived_from.append(f"adme:{value.record_id}")
        return fraction
    fraction = _explicit_float(value, label=label)
    if not 0.0 < fraction <= 1.0:
        raise AnalysisError(f"{label} must be a fraction in (0, 1], got {fraction}.")
    note = " (no binding correction)" if fraction == 1.0 else ""
    assumptions.append(f"{label}={fraction:g} supplied as an explicit scalar assumption{note}, not a measured record.")
    return fraction


def _ratio_input(
    value: InVitroAdmeRecord | float,
    *,
    compound_id: str,
    assumptions: list[str],
    derived_from: list[str],
) -> float:
    if isinstance(value, InVitroAdmeRecord):
        ratio = _measured_record_value(
            value, expected_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO, compound_id=compound_id
        )
        derived_from.append(f"adme:{value.record_id}")
        return ratio
    ratio = _explicit_float(value, label="blood_to_plasma_ratio")
    if ratio <= 0.0:
        raise AnalysisError(f"blood_to_plasma_ratio must be positive, got {ratio}.")
    note = " (blood and plasma assumed equal)" if ratio == 1.0 else ""
    assumptions.append(
        f"blood_to_plasma_ratio={ratio:g} supplied as an explicit scalar assumption{note}, not a measured record."
    )
    return ratio


def _measured_record_value(record: InVitroAdmeRecord, *, expected_kind: AdmeAssayKind, compound_id: str) -> float:
    if record.assay_kind is not expected_kind:
        raise AnalysisError(f"{record.record_id} is {record.assay_kind.value}; expected {expected_kind.value}.")
    if record.compound_id != compound_id:
        raise AnalysisError(
            f"{record.record_id} belongs to compound {record.compound_id}, not {compound_id}. "
            "Binding records are not borrowed across compounds."
        )
    canonical = record.canonical_value()
    if canonical is None or record.status is not AdmeStatus.MEASURED:
        raise AnalysisError(f"{record.record_id} is {record.status.value}; only a measured value can be used.")
    return canonical[0]


def _explicit_float(value: object, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise AnalysisError(f"{label} must be a measured InVitroAdmeRecord or an explicit number.")
    result = float(value)
    if result != result or result in {float("inf"), float("-inf")}:
        raise AnalysisError(f"{label} must be finite.")
    return result
