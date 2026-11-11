"""PharmacologyEvidence built from exposure and analysis records.

Estimated NCA parameters, PK compartment fits, direct/sigmoid Emax parameters,
in vitro ADME rows, toxicology endpoints, and exploratory IVIVE predictions
keep their evidence grade. Model effect predictions and IVIVE clearances stay
predicted. Extravascular apparent parameters (CL/F,
V/F, Vz/F) remain distinct from absolute IV clearance and volume. Nothing in
this module imports TheraWeave.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from pharmacoweave.core.evidence import (
    RESEARCH_USE_LIMITATION,
    EvidenceExchange,
    EvidenceQc,
    EvidenceRecord,
    EvidenceType,
    ProjectName,
)
from pharmacoweave.core.modeling import (
    DiagnosticStatus,
    NcaParameterSet,
    ParameterSemantic,
    auc_inf_is_reportable,
)
from pharmacoweave.dmpk.physiology import IviveClearancePrediction
from pharmacoweave.dmpk.records import AdmeStatus, InVitroAdmeRecord
from pharmacoweave.integrations.guards import screen_exchange
from pharmacoweave.pd_models.emax import DirectEmaxFit
from pharmacoweave.pk_models.fit_extravascular import ExtravascularFit
from pharmacoweave.pk_models.fit_infusion import InfusionFit
from pharmacoweave.pk_models.fit_iv import IvBolusFit
from pharmacoweave.pk_models.fit_three_compartment import ThreeCompartmentFit
from pharmacoweave.pk_models.fit_two_compartment import TwoCompartmentFit
from pharmacoweave.toxicology.records import ToxEndpoint

PkFitResult = IvBolusFit | InfusionFit | ExtravascularFit | TwoCompartmentFit | ThreeCompartmentFit

_NCA_MEASUREMENT = {
    "auc_last": "auc",
    "auc_inf": "auc",
    "auc_tau": "auc",
    "cmax": "cmax",
    "tmax": "tmax",
    "lambda_z": "lambda_z",
    "half_life": "half_life",
    "clearance": "clearance",
    "CL_over_F": "apparent_clearance",
    "volume": "volume",
    "Vz_over_F": "apparent_volume",
    "V_over_F": "apparent_volume",
}
_PK_FIT_MEASUREMENT = {
    ParameterSemantic.CL: "clearance",
    ParameterSemantic.CL_OVER_F: "apparent_clearance",
    ParameterSemantic.V: "volume",
    ParameterSemantic.V_OVER_F: "apparent_volume",
    ParameterSemantic.VZ: "volume",
    ParameterSemantic.VZ_OVER_F: "apparent_volume",
    ParameterSemantic.V1: "volume",
    ParameterSemantic.V2: "volume",
    ParameterSemantic.V3: "volume",
    ParameterSemantic.Q: "clearance",
    ParameterSemantic.Q2: "clearance",
    ParameterSemantic.Q3: "clearance",
    ParameterSemantic.KA: "ka",
    ParameterSemantic.KEL: "lambda_z",
    ParameterSemantic.K12: "k12",
    ParameterSemantic.K21: "k21",
    ParameterSemantic.HALF_LIFE: "half_life",
    ParameterSemantic.F: "bioavailability",
}
_EMAX_MEASUREMENT = {
    "Emax": "emax",
    "EC50": "ec50",
    "hill": "hill",
}


def evidence_from_nca(parameters: NcaParameterSet, *, run_id: str) -> list[EvidenceRecord]:
    """Estimated exposure parameters that are allowed to be reported.

    AUC∞ is omitted when ``auc_inf_is_reportable`` is false, even if the
    number is stored on the NCA record. Extravascular routes emit
    ``CL_over_F`` and ``Vz_over_F`` instead of ``clearance`` and ``volume``.
    """

    if not run_id.strip():
        raise ValueError("run_id is required.")
    claims: list[tuple[str, float, str]] = []
    _append(claims, "auc_last", parameters.auc_last, parameters.auc_last_unit)
    _append(claims, "auc_tau", parameters.auc_tau, parameters.auc_tau_unit)
    if auc_inf_is_reportable(parameters):
        _append(claims, "auc_inf", parameters.auc_inf, parameters.auc_inf_unit)
    _append(claims, "cmax", parameters.cmax, parameters.cmax_unit)
    _append(claims, "tmax", parameters.tmax, parameters.tmax_unit)
    _append(claims, "lambda_z", parameters.lambda_z, parameters.lambda_z_unit)
    _append(claims, "half_life", parameters.half_life, parameters.half_life_unit)
    cl_claim = (
        "CL_over_F"
        if parameters.clearance_semantic is ParameterSemantic.CL_OVER_F
        else "clearance"
    )
    vol_claim = (
        parameters.volume_semantic.value
        if parameters.volume_semantic in {ParameterSemantic.VZ_OVER_F, ParameterSemantic.V_OVER_F}
        else "volume"
    )
    _append(claims, cl_claim, parameters.clearance, parameters.clearance_unit)
    _append(claims, vol_claim, parameters.volume, parameters.volume_unit)
    derived = [
        f"dose_event:{parameters.dose_event_id}",
        f"analyte:{parameters.analyte_id}",
        f"subject:{parameters.subject_id}",
    ]
    if parameters.regimen == "single_dose":
        nca_note = (
            f"Estimated from single-dose {parameters.route.value} NCA. "
            "Not an observed concentration and not a dosing recommendation."
        )
    else:
        nca_note = (
            f"Estimated from {parameters.regimen} {parameters.route.value} NCA over one dosing interval. "
            "Steady-state clearance is Dose/AUC_tau only when the caller declared steady_state. "
            "Not an observed concentration and not a dosing recommendation."
        )
    records: list[EvidenceRecord] = []
    for claim, value, unit in claims:
        records.append(
            _record(
                evidence_id=f"pw-{run_id}-{parameters.subject_id}-{parameters.analyte_id}-{claim}",
                run_id=run_id,
                entity_type="NcaParameterSet",
                entity_id=f"{parameters.subject_id}:{parameters.analyte_id}",
                measurement_type=_NCA_MEASUREMENT[claim],
                evidence_type=EvidenceType.ESTIMATED,
                qc_status=EvidenceQc.NOT_REVIEWED,
                claim_type=claim,
                unit=unit,
                numeric_value=value,
                limitations=[RESEARCH_USE_LIMITATION, nca_note],
                derived_from=[*derived, f"nca:{parameters.subject_id}:{parameters.analyte_id}:{claim}"],
                applicability="Research-use exposure estimate.",
            )
        )
    return records


def evidence_from_pk_fit(result: PkFitResult, *, run_id: str) -> list[EvidenceRecord]:
    """Estimated compartmental PK parameters with explicit ``derived_from`` provenance."""

    if not run_id.strip():
        raise ValueError("run_id is required.")
    entity_type = type(result).__name__
    derived = [f"fit:{result.fit.fit_id}", f"subject:{result.subject_id}"]
    derived.extend(f"observation:{obs_id}" for obs_id in result.observation_ids)
    limitations = list(result.diagnostic.limitations)
    if RESEARCH_USE_LIMITATION not in limitations:
        limitations.insert(0, RESEARCH_USE_LIMITATION)
    if result.diagnostic.status is DiagnosticStatus.FAIL or not result.estimates:
        failed_measurement = (
            "apparent_clearance" if isinstance(result, ExtravascularFit) else "clearance"
        )
        return [
            _record(
                evidence_id=f"pkfit-{run_id}-{result.fit.fit_id}-not-identifiable",
                run_id=run_id,
                entity_type=entity_type,
                entity_id=result.fit.fit_id,
                measurement_type=failed_measurement,
                evidence_type=EvidenceType.ESTIMATED,
                qc_status=EvidenceQc.NOT_REVIEWED,
                claim_type="pk_fit_not_identifiable",
                unit=None,
                numeric_value=None,
                limitations=limitations,
                derived_from=derived,
                applicability="Research-use PK compartment fit. Parameters were not reported because the fit failed.",
                subject_id=result.subject_id,
            )
        ]
    return [
        _record(
            evidence_id=f"pkfit-{run_id}-{result.fit.fit_id}-{estimate.semantic.value}",
            run_id=run_id,
            entity_type=entity_type,
            entity_id=result.fit.fit_id,
            measurement_type=_PK_FIT_MEASUREMENT[estimate.semantic],
            evidence_type=EvidenceType.ESTIMATED,
            qc_status=EvidenceQc.NOT_REVIEWED,
            claim_type=estimate.semantic.value,
            unit=estimate.unit,
            numeric_value=estimate.value,
            limitations=limitations,
            derived_from=derived,
            applicability="Research-use compartmental PK parameter estimate. Not a dosing recommendation.",
            subject_id=result.subject_id,
        )
        for estimate in result.estimates
    ]


def evidence_from_urine(result: Any, *, run_id: str) -> list[EvidenceRecord]:
    """Estimated urine amount and renal clearance. Urine concentration is not exported."""

    if not run_id.strip():
        raise ValueError("run_id is required.")
    limitations = list(result.limitations)
    if RESEARCH_USE_LIMITATION not in limitations:
        limitations.insert(0, RESEARCH_USE_LIMITATION)
    limitations.extend(result.warnings)
    derived = [
        f"dose_event:{result.dose_event_id}",
        f"subject:{result.subject_id}",
        f"analyte:{result.analyte_id}",
    ]
    derived.extend(f"observation:{item.observation_id}" for item in result.intervals)
    if result.amount_excreted_mg is None and result.renal_clearance_l_h is None:
        return [
            _record(
                evidence_id=f"urine-{run_id}-{result.subject_id}-{result.analyte_id}-incomplete",
                run_id=run_id,
                entity_type="UrineExcretionResult",
                entity_id=f"{result.subject_id}:{result.analyte_id}",
                measurement_type="renal_clearance",
                evidence_type=EvidenceType.ESTIMATED,
                qc_status=EvidenceQc.NOT_REVIEWED,
                claim_type="urine_excretion_incomplete",
                unit=None,
                numeric_value=None,
                limitations=limitations,
                derived_from=derived,
                applicability="Research-use urine excretion. Amount and renal clearance were withheld.",
                subject_id=result.subject_id,
            )
        ]
    claims: list[tuple[str, float, str, str]] = []
    if result.amount_excreted_mg is not None:
        claims.append(("amount_excreted", result.amount_excreted_mg, "mg", "amount_excreted"))
    if result.renal_clearance_l_h is not None:
        claims.append(("renal_clearance", result.renal_clearance_l_h, "L/h", "renal_clearance"))
    return [
        _record(
            evidence_id=f"urine-{run_id}-{result.subject_id}-{result.analyte_id}-{claim}",
            run_id=run_id,
            entity_type="UrineExcretionResult",
            entity_id=f"{result.subject_id}:{result.analyte_id}",
            measurement_type=measurement,
            evidence_type=EvidenceType.ESTIMATED,
            qc_status=EvidenceQc.NOT_REVIEWED,
            claim_type=claim,
            unit=unit,
            numeric_value=value,
            limitations=limitations,
            derived_from=derived,
            applicability=(
                "Research-use urine excretion estimate. "
                "Renal clearance uses a caller-supplied plasma AUC and is not plasma clearance."
            ),
            subject_id=result.subject_id,
        )
        for claim, value, unit, measurement in claims
    ]


def evidence_from_emax(result: DirectEmaxFit, *, run_id: str) -> list[EvidenceRecord]:
    """Estimated Emax parameters, plus predicted effects when the fit reports them."""

    if not run_id.strip():
        raise ValueError("run_id is required.")
    derived = [f"fit:{result.fit.fit_id}", *[f"effect:{item}" for item in result.effect_ids]]
    derived.extend(f"exposure:{item}" for item in result.exposure_refs)
    if result.diagnostic.status.value == "fail" or not result.estimates:
        return [
            _record(
                evidence_id=f"{result.fit.fit_id}-not-identifiable",
                run_id=run_id,
                entity_type="DirectEmaxFit",
                entity_id=result.fit.fit_id,
                measurement_type="emax",
                evidence_type=EvidenceType.ESTIMATED,
                qc_status=EvidenceQc.NOT_REVIEWED,
                claim_type="direct_emax_not_identifiable",
                unit=None,
                numeric_value=None,
                limitations=list(result.diagnostic.limitations),
                derived_from=derived,
                applicability="Research-use direct Emax. EC50 was not reported because the fit failed.",
            )
        ]
    hill_fixed = any(est.name == "hill" and est.fixed for est in result.estimates)
    param_note = (
        "E0 was fixed and hill was fixed at 1. This is not an observed effect."
        if hill_fixed
        else "E0 was fixed and hill was estimated on the original effect scale. This is not an observed effect."
    )
    records = [
        _record(
            evidence_id=f"{result.fit.fit_id}-{estimate.name}",
            run_id=run_id,
            entity_type="DirectEmaxFit",
            entity_id=result.fit.fit_id,
            measurement_type=_EMAX_MEASUREMENT.get(estimate.name, estimate.name.lower()),
            evidence_type=EvidenceType.ESTIMATED,
            qc_status=EvidenceQc.NOT_REVIEWED,
            claim_type=estimate.name,
            unit=estimate.unit,
            numeric_value=estimate.value,
            limitations=[
                RESEARCH_USE_LIMITATION,
                param_note,
            ],
            derived_from=derived,
            applicability="Research-use direct-effect parameter. Not a dosing recommendation.",
        )
        for estimate in result.estimates
        if not estimate.fixed
    ]
    records.extend(
        _record(
            evidence_id=f"{result.fit.fit_id}-pred-{point.effect_id}",
            run_id=run_id,
            entity_type="PredictedEffect",
            entity_id=point.effect_id,
            measurement_type="effect",
            evidence_type=EvidenceType.PREDICTED,
            qc_status=EvidenceQc.NOT_APPLICABLE,
            claim_type="direct_emax_prediction",
            unit=point.unit,
            numeric_value=point.value,
            limitations=[
                RESEARCH_USE_LIMITATION,
                "Model prediction from direct Emax. Not an observed effect.",
            ],
            derived_from=list(point.derived_from),
            applicability="Research-use model prediction.",
        )
        for point in result.predicted
    )
    return records


def evidence_from_adme(record: InVitroAdmeRecord, *, run_id: str) -> EvidenceRecord | None:
    """Observed in vitro evidence. A not-measured row does not become a number.

    ``measurement_type`` is the assay kind, not ``clearance``, so a microsomal
    result is not labeled as a plasma clearance.
    """

    if record.status is AdmeStatus.NOT_MEASURED:
        return None
    if not run_id.strip():
        raise ValueError("run_id is required.")
    failed = record.status is AdmeStatus.FAILED_QC
    limitations = list(record.limitations)
    if failed and record.qc_reason and record.qc_reason not in limitations:
        limitations.append(record.qc_reason)
    return _record(
        evidence_id=f"adme-{run_id}-{record.record_id}",
        run_id=run_id,
        entity_type="InVitroAdmeRecord",
        entity_id=record.record_id,
        measurement_type=record.assay_kind.value,
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.REJECTED if failed else EvidenceQc.NOT_REVIEWED,
        claim_type=f"{record.assay_kind.value}_failed_qc" if failed else record.assay_kind.value,
        unit=None if failed else record.unit,
        numeric_value=None if failed else record.value,
        limitations=limitations,
        derived_from=[f"adme:{record.record_id}", f"compound:{record.compound_id}"],
        applicability="In vitro ADME observation. Not an in vivo clearance and not a predicted ADMET value.",
    )


def evidence_from_tox(record: ToxEndpoint, *, run_id: str) -> EvidenceRecord:
    """Observed toxicology evidence. The applicability text refuses clinical dosing advice."""

    if not run_id.strip():
        raise ValueError("run_id is required.")
    derived = [f"tox:{record.endpoint_id}", f"compound:{record.compound_id}"]
    if record.dose_event_id is not None:
        derived.append(f"dose_event:{record.dose_event_id}")
    if record.exposure_ref:
        derived.append(f"exposure:{record.exposure_ref}")
    return _record(
        evidence_id=f"tox-{run_id}-{record.endpoint_id}",
        run_id=run_id,
        entity_type="ToxEndpoint",
        entity_id=record.endpoint_id,
        measurement_type="toxicity_endpoint",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.NOT_REVIEWED,
        claim_type=record.kind.value,
        unit=record.unit,
        numeric_value=record.numeric_value,
        limitations=list(record.limitations),
        derived_from=derived,
        applicability="Research-use toxicology observation. Not a clinical dosing recommendation.",
        subject_id=record.subject_id,
    )


def evidence_from_ivive(prediction: IviveClearancePrediction, *, run_id: str) -> EvidenceRecord:
    """Predicted hepatic plasma clearance from an exploratory IVIVE.

    The record is ``predicted`` with ``qc_status=not_applicable``. It cites
    the CLint, binding, and physiology inputs in ``derived_from`` and is never
    labeled observed or estimated.
    """

    if not run_id.strip():
        raise ValueError("run_id is required.")
    limitations = list(prediction.limitations)
    limitations.extend(f"Assumption: {item}" for item in prediction.assumptions)
    return _record(
        evidence_id=f"ivive-{run_id}-{prediction.prediction_id}",
        run_id=run_id,
        entity_type="CandidateCompound",
        entity_id=prediction.compound_id,
        measurement_type="clearance",
        evidence_type=EvidenceType.PREDICTED,
        qc_status=EvidenceQc.NOT_APPLICABLE,
        claim_type="ivive_hepatic_plasma_clearance",
        unit=prediction.clearance_unit,
        numeric_value=prediction.hepatic_plasma_clearance_l_h,
        limitations=limitations,
        derived_from=[f"ivive:{prediction.prediction_id}", *prediction.derived_from],
        applicability=(
            f"Exploratory {prediction.model.value} hepatic IVIVE for {prediction.species}. "
            "Prior for calibration audit only; not a dosing recommendation."
        ),
    )


def collect_run_evidence(
    *,
    run_id: str,
    nca_parameters: Sequence[NcaParameterSet] | None = None,
    pk_fits: Sequence[PkFitResult] | None = None,
    emax_fits: Sequence[DirectEmaxFit] | None = None,
    adme_records: Sequence[InVitroAdmeRecord] | None = None,
    tox_endpoints: Sequence[ToxEndpoint] | None = None,
    ivive_predictions: Sequence[IviveClearancePrediction] | None = None,
    urine_results: Sequence[Any] | None = None,
) -> list[EvidenceRecord]:
    """Build and screen PharmacologyEvidence records across all run outputs."""

    records: list[EvidenceRecord] = []
    for nca in nca_parameters or ():
        records.extend(evidence_from_nca(nca, run_id=run_id))
    for pk_fit in pk_fits or ():
        records.extend(evidence_from_pk_fit(pk_fit, run_id=run_id))
    for emax_fit in emax_fits or ():
        records.extend(evidence_from_emax(emax_fit, run_id=run_id))
    for adme in adme_records or ():
        built = evidence_from_adme(adme, run_id=run_id)
        if built is not None:
            records.append(built)
    for tox in tox_endpoints or ():
        records.append(evidence_from_tox(tox, run_id=run_id))
    for prediction in ivive_predictions or ():
        records.append(evidence_from_ivive(prediction, run_id=run_id))
    for urine in urine_results or ():
        records.extend(evidence_from_urine(urine, run_id=run_id))
    if records:
        screen_exchange(EvidenceExchange(records=records))
    return records


def exchange_from_records(records: list[EvidenceRecord]) -> EvidenceExchange:
    """Screen and return an exchange. Grades are not rewritten."""

    exchange = EvidenceExchange(records=records)
    screen_exchange(exchange)
    return exchange


def round_trip_exchange(exchange: EvidenceExchange) -> EvidenceExchange:
    """Dump an exchange to JSON and load it again, then screen it."""

    loaded = EvidenceExchange.model_validate_json(exchange.model_dump_json())
    screen_exchange(loaded)
    return loaded


def write_exchange(exchange: EvidenceExchange, path: str | Path) -> None:
    """Write one screened evidence exchange. The file does not require TheraWeave."""

    screen_exchange(exchange)
    destination = Path(path)
    destination.write_text(exchange.model_dump_json(indent=2) + "\n", encoding="utf-8")


def _append(claims: list[tuple[str, float, str]], name: str, value: float | None, unit: str | None) -> None:
    if value is None or unit is None:
        return
    claims.append((name, value, unit))


def _record(
    *,
    evidence_id: str,
    run_id: str,
    entity_type: str,
    entity_id: str,
    measurement_type: str,
    evidence_type: EvidenceType,
    qc_status: EvidenceQc,
    claim_type: str,
    unit: str | None,
    numeric_value: float | None,
    limitations: list[str],
    derived_from: list[str],
    applicability: str,
    subject_id: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        evidence_id=evidence_id,
        source_project=ProjectName.PHARMACOWEAVE,
        source_run_id=run_id,
        entity_type=entity_type,
        entity_id=entity_id,
        measurement_type=measurement_type,
        evidence_type=evidence_type,
        qc_status=qc_status,
        claim_type=claim_type,
        unit=unit,
        numeric_value=numeric_value,
        limitations=limitations,
        derived_from=derived_from,
        applicability=applicability,
        subject_id=subject_id,
    )
