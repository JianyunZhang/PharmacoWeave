"""Structural model specifications that have a tested equation."""

from pharmacoweave.core.modeling import ErrorModel, ModelSpec

ONE_COMPARTMENT_IV_BOLUS = ModelSpec(
    model_id="one_compartment_iv_bolus",
    model_type="one_compartment_iv_bolus",
    equations_ref="C(t) = (Dose/V) * exp(-(CL/V)*t)",
    parameters=["CL", "V"],
    error_model=ErrorModel.PROPORTIONAL,
    assumptions=[
        "Input is an instantaneous IV bolus.",
        "Elimination is linear and first order.",
        "The observation is the central-compartment concentration.",
        "CL and V can be estimated by nonlinear least squares for this structure.",
        "The default objective uses proportional residuals. Additive residuals are optional.",
        "Non-quantified samples are omitted from the objective.",
    ],
    question_of_interest="What CL and V does a one-compartment IV bolus assign to these concentrations?",
    context_of_use=(
        "Individual fit and reference equation for single-dose IV bolus concentrations. "
        "Research use only. Not a dosing recommendation."
    ),
    estimation_implemented=True,
)

BUILTIN_SPECS = {"one_compartment_iv_bolus": ONE_COMPARTMENT_IV_BOLUS}
