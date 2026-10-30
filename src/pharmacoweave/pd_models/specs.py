"""Direct-effect model specifications that have a tested equation."""

from pharmacoweave.core.modeling import ErrorModel, ModelSpec

DIRECT_EMAX = ModelSpec(
    model_id="direct_emax_hill1",
    model_type="direct_emax",
    equations_ref="E = E0 + Emax * C / (EC50 + C)",
    parameters=["E0", "Emax", "EC50", "hill"],
    error_model=ErrorModel.UNSPECIFIED,
    assumptions=[
        "The baseline E0 is fixed and is not estimated.",
        "The Hill coefficient is fixed at 1.",
        "Each effect is paired with an explicit positive concentration.",
        "AUC and Cmax are not substituted for that concentration.",
        "No lag and no indirect response are estimated.",
    ],
    question_of_interest=(
        "What Emax and EC50 does a direct concentration-effect curve support when E0 and hill are fixed?"
    ),
    context_of_use=(
        "Research-use direct-effect Emax for a simple concentration-effect curve. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)
