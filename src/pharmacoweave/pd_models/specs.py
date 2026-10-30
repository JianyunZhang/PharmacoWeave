"""Direct-effect and sigmoid Emax model specifications with tested equations."""

from pharmacoweave.core.modeling import ErrorModel, ModelSpec

DIRECT_EMAX = ModelSpec(
    model_id="direct_emax_hill1",
    model_type="direct_emax",
    equations_ref="E = E0 + Emax * C / (EC50 + C)",
    parameters=["E0", "Emax", "EC50", "hill"],
    error_model=ErrorModel.ADDITIVE,
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

SIGMOID_EMAX = ModelSpec(
    model_id="sigmoid_emax",
    model_type="sigmoid_emax",
    equations_ref="E = E0 + Emax * C^hill / (EC50^hill + C^hill)",
    parameters=["E0", "Emax", "EC50", "hill"],
    error_model=ErrorModel.ADDITIVE,
    assumptions=[
        "The baseline E0 is fixed and is not estimated.",
        "The Hill coefficient is positive and requires observations across the transition zone.",
        "Each effect is paired with an explicit positive concentration.",
        "AUC and Cmax are not substituted for that concentration.",
        "No lag and no indirect response are estimated.",
    ],
    question_of_interest=(
        "What Emax, EC50, and Hill coefficient does a direct concentration-effect curve support when E0 is fixed?"
    ),
    context_of_use=(
        "Research-use sigmoid direct-effect Emax for concentration-response curves with transition support. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)
