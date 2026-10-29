"""Record-level scientific rules."""

from datetime import UTC, datetime

import pytest
from pydantic import ValidationError

from pharmacoweave.core.entities import (
    AnalyteForm,
    DoseEvent,
    Observation,
    ObservationStatus,
    ObservationType,
    Route,
)
from pharmacoweave.core.modeling import (
    DiagnosticRecord,
    DiagnosticStatus,
    ModelFit,
    NcaParameterSet,
    ParameterEstimate,
    ParameterSemantic,
    auc_inf_is_reportable,
)


def test_blq_cannot_store_a_number() -> None:
    with pytest.raises(ValidationError, match="numeric_value"):
        Observation(
            observation_id="OBS-1",
            sample_id="SAMPLE-1",
            analyte_id="ANALYTE-1",
            observation_type=ObservationType.CONCENTRATION,
            status=ObservationStatus.BLQ,
            numeric_value=0.0,
            unit="ug/mL",
            assay_id="ASSAY-1",
            analyte_form=AnalyteForm.PARENT,
        )


def test_quantified_concentration_cannot_be_zero() -> None:
    with pytest.raises(ValidationError, match="positive"):
        Observation(
            observation_id="OBS-1",
            sample_id="SAMPLE-1",
            analyte_id="ANALYTE-1",
            observation_type=ObservationType.CONCENTRATION,
            status=ObservationStatus.QUANTIFIED,
            numeric_value=0.0,
            unit="ug/mL",
            assay_id="ASSAY-1",
        )


def test_naive_dose_time_is_rejected() -> None:
    with pytest.raises(ValidationError, match="timezone"):
        DoseEvent(
            dose_event_id="DOSE-1",
            subject_id="SUBJ-1",
            compound_id="CMPD-1",
            absolute_time=datetime(2024, 1, 15, 9, 0),
            dose_value=1.0,
            dose_unit="mg",
            route=Route.IV_BOLUS,
        )


def test_infusion_requires_duration() -> None:
    with pytest.raises(ValidationError, match="infusion"):
        DoseEvent(
            dose_event_id="DOSE-1",
            subject_id="SUBJ-1",
            compound_id="CMPD-1",
            absolute_time=datetime(2024, 1, 15, 9, 0, tzinfo=UTC),
            dose_value=1.0,
            dose_unit="mg",
            route=Route.IV_INFUSION,
        )


def test_intravenous_clearance_is_not_apparent() -> None:
    with pytest.raises(ValidationError, match="CL_over_F"):
        ParameterEstimate(
            fit_id="FIT-1",
            parameter="CL",
            semantic=ParameterSemantic.CL_OVER_F,
            value=0.05,
            unit="L/h",
            estimation_method="reference",
            route=Route.IV_BOLUS,
        )


def test_oral_clearance_must_be_apparent() -> None:
    with pytest.raises(ValidationError, match="CL/F"):
        ParameterEstimate(
            fit_id="FIT-1",
            parameter="CL",
            semantic=ParameterSemantic.CL,
            value=0.05,
            unit="L/h",
            estimation_method="reference",
            route=Route.ORAL,
        )
    apparent = ParameterEstimate(
        fit_id="FIT-1",
        parameter="CL/F",
        semantic=ParameterSemantic.CL_OVER_F,
        value=0.05,
        unit="L/h",
        estimation_method="reference",
        route=Route.ORAL,
    )
    assert apparent.semantic is ParameterSemantic.CL_OVER_F


def test_unspecified_route_cannot_carry_clearance() -> None:
    with pytest.raises(ValidationError, match="Route other"):
        ParameterEstimate(
            fit_id="FIT-1",
            parameter="CL",
            semantic=ParameterSemantic.CL,
            value=0.05,
            unit="L/h",
            estimation_method="reference",
            route=Route.OTHER,
        )


def test_terminal_parameters_require_a_qualified_phase() -> None:
    with pytest.raises(ValidationError, match="terminal phase"):
        NcaParameterSet(
            subject_id="SUBJ-1",
            analyte_id="ANALYTE-1",
            route=Route.IV_BOLUS,
            dose_event_id="DOSE-1",
            auc_method="linear_up_log_down",
            terminal_phase_qualified=False,
            auc_inf=5.0,
            auc_inf_unit="mg*h/L",
            auc_inf_extrap_fraction=0.1,
        )


def test_auc_infinity_reporting_flag() -> None:
    record = NcaParameterSet(
        subject_id="SUBJ-1",
        analyte_id="ANALYTE-1",
        route=Route.ORAL,
        dose_event_id="DOSE-1",
        auc_method="linear_up_log_down",
        terminal_phase_qualified=True,
        auc_inf=5.0,
        auc_inf_unit="mg*h/L",
        auc_inf_extrap_fraction=0.25,
        clearance_semantic=ParameterSemantic.CL_OVER_F,
        clearance=1.0,
        clearance_unit="L/h",
    )
    assert auc_inf_is_reportable(record, max_extrapolated_fraction=0.20) is False
    assert auc_inf_is_reportable(record, max_extrapolated_fraction=0.30) is True


def test_passing_diagnostic_states_limitations() -> None:
    with pytest.raises(ValidationError, match="limitations"):
        DiagnosticRecord(fit_id="FIT-1", status=DiagnosticStatus.PASS, limitations=[])
    with pytest.raises(ValidationError, match="converged"):
        ModelFit(
            fit_id="FIT-1",
            model_id="one_compartment_iv_bolus",
            dataset_id="DATA-1",
            software_version="0.0.1",
            estimation_method="not-run",
            context_of_use="Reference fit record for the validation test.",
            diagnostic_status=DiagnosticStatus.PASS,
            converged=False,
            limitations=["The fit did not converge."],
        )
