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
        "CL and V may be estimated by log-linear ordinary least squares or by weighted nonlinear least squares.",
        "Supported weights are 1/y and 1/y^2. BLQ observations are omitted and are not replaced with zero.",
        "Two-compartment, infusion, and extravascular models are not estimated.",
    ],
    question_of_interest="What CL and V does a one-compartment IV bolus concentration curve support?",
    context_of_use=(
        "Research-use estimation of CL and V for a single-dose one-compartment IV bolus. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)

BUILTIN_SPECS = {"one_compartment_iv_bolus": ONE_COMPARTMENT_IV_BOLUS}
