"""Structural model specifications that have a tested equation."""

from pharmacoweave.core.modeling import ErrorModel, ModelSpec

ONE_COMPARTMENT_IV_BOLUS = ModelSpec(
    model_id="one_compartment_iv_bolus",
    model_type="one_compartment_iv_bolus",
    equations_ref="C(t) = (Dose/V) * exp(-(CL/V)*t)",
    parameters=["CL", "V"],
    error_model=ErrorModel.UNSPECIFIED,
    assumptions=[
        "Input is an instantaneous IV bolus.",
        "Elimination is linear and first order.",
        "The observation is the central-compartment concentration.",
        "This specification does not estimate parameters.",
    ],
    question_of_interest="What concentration does a known one-compartment IV bolus produce?",
    context_of_use=(
        "Reference equation for single-dose IV bolus concentrations. "
        "Not a fitted model and not a dosing recommendation."
    ),
    estimation_implemented=False,
)

BUILTIN_SPECS = {"one_compartment_iv_bolus": ONE_COMPARTMENT_IV_BOLUS}
