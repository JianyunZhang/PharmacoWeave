"""Structural PK model specifications that have a tested equation and estimator."""

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
        "Supported weights are uniform, 1/y, and 1/y^2. BLQ observations are omitted and are not replaced with zero.",
    ],
    question_of_interest="What CL and V does a one-compartment IV bolus concentration curve support?",
    context_of_use=(
        "Research-use estimation of CL and V for a single-dose one-compartment IV bolus. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)

ONE_COMPARTMENT_IV_INFUSION = ModelSpec(
    model_id="one_compartment_iv_infusion",
    model_type="one_compartment_iv_infusion",
    equations_ref=(
        "C(t) = (Rate/CL) * (1 - exp(-(CL/V)*t)) during a constant-rate IV infusion, "
        "then C(T) * exp(-(CL/V)*(t-T)); Rate = Dose/T"
    ),
    parameters=["CL", "V"],
    error_model=ErrorModel.UNSPECIFIED,
    assumptions=[
        "Input is a single-dose constant-rate IV infusion of a declared positive duration.",
        "The infusion duration is shorter than the last quantified time.",
        "Elimination is linear and first order.",
        "The observation is the central-compartment concentration.",
        "CL and V are absolute intravenous parameters, not CL/F or V/F.",
        "Supported weights are uniform, 1/y, and 1/y^2. BLQ observations are omitted and are not replaced with zero.",
    ],
    question_of_interest="What CL and V does a one-compartment constant-rate IV infusion curve support?",
    context_of_use=(
        "Research-use estimation of CL and V for a single-dose one-compartment constant-rate IV infusion. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)

ONE_COMPARTMENT_EXTRAVASCULAR = ModelSpec(
    model_id="one_compartment_extravascular",
    model_type="one_compartment_extravascular",
    equations_ref="C(t) = (Dose/(V/F)) * (ka/(ka - (CL/F)/(V/F))) * (exp(-((CL/F)/(V/F))*t) - exp(-ka*t))",
    parameters=["ka", "CL/F", "V/F"],
    error_model=ErrorModel.UNSPECIFIED,
    assumptions=[
        "Input is single-dose extravascular first-order absorption without lag time.",
        "Absorption rate exceeds elimination rate (ka > ke); flip-flop kinetics cannot be separated without an IV arm.",
        "Bioavailability F is unknown in a single-arm extravascular design, so apparent CL/F and V/F are reported.",
        "Supported weights are uniform, 1/y, and 1/y^2. BLQ observations are omitted and are not replaced with zero.",
    ],
    question_of_interest="What ka, CL/F, and V/F does a single-dose extravascular concentration curve support?",
    context_of_use=(
        "Research-use estimation of ka, CL/F, and V/F for a single-dose first-order absorption curve. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)

TWO_COMPARTMENT_IV_BOLUS = ModelSpec(
    model_id="two_compartment_iv_bolus",
    model_type="two_compartment_iv_bolus",
    equations_ref="C(t) = A*exp(-alpha*t) + B*exp(-beta*t) parameterized by (CL, V1, Q, V2)",
    parameters=["CL", "V1", "Q", "V2"],
    error_model=ErrorModel.UNSPECIFIED,
    assumptions=[
        "Input is an instantaneous IV bolus into the central compartment.",
        "Distribution between central (V1) and peripheral (V2) compartments is linear and governed by Q.",
        "Elimination occurs linearly from the central compartment (CL).",
        "Supported weights are uniform, 1/y, and 1/y^2. BLQ observations are omitted and are not replaced with zero.",
    ],
    question_of_interest="What CL, V1, Q, and V2 does a biexponential IV bolus concentration curve support?",
    context_of_use=(
        "Research-use estimation of CL, V1, Q, and V2 for a single-dose two-compartment IV bolus. "
        "Not a dosing recommendation and not a regulatory validation."
    ),
    estimation_implemented=True,
)

BUILTIN_SPECS = {
    "one_compartment_iv_bolus": ONE_COMPARTMENT_IV_BOLUS,
    "one_compartment_iv_infusion": ONE_COMPARTMENT_IV_INFUSION,
    "one_compartment_extravascular": ONE_COMPARTMENT_EXTRAVASCULAR,
    "two_compartment_iv_bolus": TWO_COMPARTMENT_IV_BOLUS,
}
