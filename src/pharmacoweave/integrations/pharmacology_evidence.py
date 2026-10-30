"""PharmacologyEvidence built from exposure and analysis records.

Estimated NCA parameters, direct-Emax parameters, in vitro ADME rows, and
toxicology endpoints keep their evidence grade. Model effect predictions stay
predicted. Nothing in this module imports TheraWeave.
"""

from pathlib import Path

from pharmacoweave.core.evidence import (
    RESEARCH_USE_LIMITATION,
    EvidenceExchange,
    EvidenceQc,
    EvidenceRecord,
    EvidenceType,
    ProjectName,
)
from pharmacoweave.core.modeling import NcaParameterSet, auc_inf_is_reportable
from pharmacoweave.dmpk.records import AdmeStatus, InVitroAdmeRecord
from pharmacoweave.integrations.guards import screen_exchange
from pharmacoweave.pd_models.emax import DirectEmaxFit
from pharmacoweave.toxicology.records import ToxEndpoint

_NCA_MEASUREMENT = {
    "auc_last": "auc",
    "auc_inf": "auc",
    "cmax": "cmax",
    "tmax": "tmax",
    "lambda_z": "lambda_z",
    "half_life": "half_life",
    "clearance": "clearance",
    "volume": "volume",
}
_NCA_NOTE = (
    "Estimated from single-dose IV bolus NCA. Not an observed concentration and not a dosing recommendation."
)


def evidence_from_nca(parameters: NcaParameterSet, *, run_id: str) -> list[EvidenceRecord]:
    """Estimated exposure parameters that are allowed to be reported.

    AUC∞ is omitted when ``auc_inf_is_reportable`` is false, even if the
    number is stored on the NCA record.
    """

    if not run_id.strip():
        raise ValueError("run_id is required.")
    claims: list[tuple[str, float, str]] = []
    _append(claims, "auc_last", parameters.auc_last, parameters.auc_last_unit)
    if auc_inf_is_reportable(parameters):
        _append(claims, "auc_inf", parameters.auc_inf, parameters.auc_inf_unit)
    _append(claims, "cmax", parameters.cmax, parameters.cmax_unit)
    _append(claims, "tmax", parameters.tmax, parameters.tmax_unit)
    _append(claims, "lambda_z", parameters.lambda_z, parameters.lambda_z_unit)
    _append(claims, "half_life", parameters.half_life, parameters.half_life_unit)
    _append(claims, "clearance", parameters.clearance, parameters.clearance_unit)
    _append(claims, "volume", parameters.volume, parameters.volume_unit)
    derived = [
        f"dose_event:{parameters.dose_event_id}",
        f"analyte:{parameters.analyte_id}",
        f"subject:{parameters.subject_id}",
    ]
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
                limitations=[RESEARCH_USE_LIMITATION, _NCA_NOTE],
                derived_from=[*derived, f"nca:{parameters.subject_id}:{parameters.analyte_id}:{claim}"],
                applicability="Research-use exposure estimate.",
            )
        )
    return records


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
    records = [
        _record(
            evidence_id=f"{result.fit.fit_id}-{estimate.name}",
            run_id=run_id,
            entity_type="DirectEmaxFit",
            entity_id=result.fit.fit_id,
            measurement_type="emax" if estimate.name == "Emax" else "ec50",
            evidence_type=EvidenceType.ESTIMATED,
            qc_status=EvidenceQc.NOT_REVIEWED,
            claim_type=estimate.name,
            unit=estimate.unit,
            numeric_value=estimate.value,
            limitations=[
                RESEARCH_USE_LIMITATION,
                "E0 was fixed and hill was fixed at 1. This is not an observed effect.",
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
    )


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
    )
