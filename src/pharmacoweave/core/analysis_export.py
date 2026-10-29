"""Analysis artifacts handed to the result-bundle writer.

The writer records these artifacts. It does not calculate NCA or fit a model.
"""

from dataclasses import dataclass, field

from pharmacoweave.core.modeling import AnalysisRecords


@dataclass(frozen=True)
class PredictionRow:
    """One observed-versus-predicted concentration used in an accepted fit."""

    fit_id: str
    observation_id: str
    subject_id: str
    time_h: float
    observed: float
    predicted: float
    residual: float
    weighted_residual: float
    unit: str


@dataclass(frozen=True)
class AnalysisPayload:
    """What an analysis run asks the bundle writer to store."""

    records: AnalysisRecords
    predictions: tuple[PredictionRow, ...] = ()
    concentration_time_svg: str = ""
    goodness_of_fit_svg: str | None = None
    residual_diagnostics_svg: str | None = None
    analysis_parameters: dict[str, object] = field(default_factory=dict)
