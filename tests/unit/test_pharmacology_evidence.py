"""PharmacologyEvidence keeps grades and derived_from links."""

import ast
from pathlib import Path

import pytest
from examples.direct_effect_and_evidence import E0, build_adme, build_effects, build_exposures, build_tox
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import Route
from pharmacoweave.core.evidence import EvidenceQc, EvidenceRecord, EvidenceType, ProjectName
from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.errors import ContractViolation
from pharmacoweave.integrations import (
    evidence_from_adme,
    evidence_from_emax,
    evidence_from_nca,
    evidence_from_tox,
    exchange_from_records,
    round_trip_exchange,
    write_exchange,
)
from pharmacoweave.integrations.guards import load_exchange, promote_to_observation, screen_record
from pharmacoweave.integrations.theraweave import attach_candidate_bundle
from pharmacoweave.pd_models import DirectEffectPair, fit_direct_emax
from pharmacoweave.pkpd import fit_linked_direct_emax


def test_nca_evidence_cites_the_dose_and_skips_unreportable_auc_inf(repo_root) -> None:
    study = build_study()
    parameters = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    records = evidence_from_nca(parameters[0], run_id="run-001")
    by_claim = {record.claim_type: record for record in records}
    assert by_claim["clearance"].evidence_type is EvidenceType.ESTIMATED
    assert by_claim["auc_inf"].evidence_type is EvidenceType.ESTIMATED
    assert "dose_event:DOSE-SUBJ-01" in by_claim["clearance"].derived_from
    assert by_claim["clearance"].source_project is ProjectName.PHARMACOWEAVE
    withheld = NcaParameterSet(
        subject_id="SUBJ-01",
        analyte_id="PW-0001-PARENT",
        route=Route.IV_BOLUS,
        dose_event_id="DOSE-SUBJ-01",
        auc_method="linear_up_log_down",
        terminal_phase_qualified=True,
        auc_last=4.0,
        auc_last_unit="mg*h/L",
        auc_inf=10.0,
        auc_inf_unit="mg*h/L",
        auc_inf_extrap_fraction=0.5,
        cmax=1.0,
        cmax_unit="mg/L",
    )
    withheld_claims = {record.claim_type for record in evidence_from_nca(withheld, run_id="run-001")}
    assert "auc_last" in withheld_claims
    assert "auc_inf" not in withheld_claims


def test_round_trip_keeps_predicted_admet_apart_from_estimates(tmp_path, repo_root) -> None:
    study = build_study()
    nca = evidence_from_nca(study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")[0], run_id="run-001")
    _linked, fit = fit_linked_direct_emax(build_effects(), build_exposures(), e0=E0, e0_unit="percent")
    emax_records = evidence_from_emax(fit, run_id="run-001")
    adme_records = [item for item in (evidence_from_adme(row, run_id="run-001") for row in build_adme()) if item]
    tox = evidence_from_tox(build_tox(), run_id="run-001")
    predicted = load_exchange(repo_root / "examples" / "cross_project_demo" / "theraweave_predicted_admet.json")
    exchange = exchange_from_records([*nca, *emax_records, *adme_records, tox, *predicted.records])
    restored = round_trip_exchange(exchange)
    grades = {(record.evidence_id, record.evidence_type, record.claim_type) for record in restored.records}
    assert ("TW-PRED-001", EvidenceType.PREDICTED, "predicted_admet") in grades
    assert any(item[1] is EvidenceType.ESTIMATED and item[2] == "clearance" for item in grades)
    assert any(item[1] is EvidenceType.PREDICTED and item[2] == "direct_emax_prediction" for item in grades)
    assert any(record.measurement_type == "microsomal_clearance" for record in adme_records)
    assert all(record.measurement_type != "clearance" for record in adme_records)
    assert evidence_from_adme(build_adme()[-1], run_id="run-001") is None
    assert "clinical dosing recommendation" in tox.applicability
    path = tmp_path / "evidence.json"
    write_exchange(restored, path)
    attached = attach_candidate_bundle(study, path)
    assert len(attached.observations) == len(study.observations)
    assert any(record.evidence_type is EvidenceType.PREDICTED for record in attached.external_evidence)
    assert any(record.evidence_type is EvidenceType.ESTIMATED for record in attached.external_evidence)
    with pytest.raises(ContractViolation):
        promote_to_observation(nca[0])


def test_failed_emax_evidence_does_not_report_ec50() -> None:
    pairs = [
        DirectEffectPair(
            effect_id=f"EFF-{index}",
            exposure_ref=f"CONC-{index}",
            concentration=float(index),
            concentration_unit="mg/L",
            effect=20.0,
            effect_unit="percent",
        )
        for index in range(1, 4)
    ]
    fit = fit_direct_emax(pairs, e0=10.0, e0_unit="percent")
    records = evidence_from_emax(fit, run_id="run-001")
    assert len(records) == 1
    assert records[0].claim_type == "direct_emax_not_identifiable"
    assert records[0].numeric_value is None
    assert records[0].evidence_type is EvidenceType.ESTIMATED
    assert records[0].derived_from
    screened = exchange_from_records(records)
    assert screened.records[0].evidence_id == records[0].evidence_id


def test_pharmacoweave_evidence_requires_derived_from_and_stays_out_of_the_observation_table() -> None:
    bare = EvidenceRecord(
        evidence_id="PW-BARE",
        source_project=ProjectName.PHARMACOWEAVE,
        source_run_id="run-001",
        entity_type="NcaParameterSet",
        entity_id="SUBJ-01",
        measurement_type="clearance",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.NOT_REVIEWED,
        claim_type="clearance",
        unit="L/h",
        numeric_value=0.05,
    )
    with pytest.raises(ContractViolation, match="derived_from"):
        screen_record(bare)
    observed = bare.model_copy(
        update={"evidence_type": EvidenceType.OBSERVED, "derived_from": ["dose_event:DOSE-1"]}
    )
    with pytest.raises(ContractViolation, match="observation table"):
        screen_record(observed)


def test_phase2_modules_do_not_import_external_engines() -> None:
    forbidden = {"pknca", "nlmixr2", "nonmem", "omicsweave", "theraweave", "ephysweave"}
    root = Path(__file__).resolve().parents[2] / "src" / "pharmacoweave"
    modules = ("pd_models", "pkpd", "dmpk", "toxicology", "integrations")
    for name in modules:
        for path in (root / name).rglob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            imported: set[str] = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0].lower() for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0].lower())
            assert imported.isdisjoint(forbidden)
