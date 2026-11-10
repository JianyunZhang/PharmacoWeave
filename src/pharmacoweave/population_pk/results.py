"""Read a caller-supplied population result table. This function starts no subprocess."""

from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Identifier, Record
from pharmacoweave.errors import AnalysisError

ParameterName = Literal["CL", "V", "CL/F", "V/F"]
EvidenceGrade = Literal["estimated", "predicted"]
_REQUIRED = ("subject_id", "parameter_name", "estimate", "unit", "estimation_method", "evidence_grade")
_FORBIDDEN = {"dosing_recommendation", "constitutes_dosing_recommendation"}


class PopulationResultRow(Record):
    """One row kept at the evidence grade the caller declared."""

    subject_id: Identifier
    parameter_name: ParameterName
    estimate: float
    unit: str = Field(min_length=1)
    estimation_method: str = Field(min_length=1)
    evidence_grade: EvidenceGrade


class PopulationResultTable(Record):
    """A checked table. Reading it does not launch an external engine."""

    rows: list[PopulationResultRow]
    limitations: list[str]
    schema_version: Literal["0.1.0"] = "0.1.0"


def load_population_result_table(rows: list[dict[str, object]]) -> PopulationResultTable:
    """Check a finished population table. Missing columns and relabeling reject the whole table."""

    if not rows:
        raise AnalysisError("A population result table requires at least one row.")
    checked: list[PopulationResultRow] = []
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            raise AnalysisError(f"Population result row {index} must be a mapping.")
        if _FORBIDDEN.intersection(row):
            raise AnalysisError("A population result table cannot carry a dosing recommendation.")
        missing = [name for name in _REQUIRED if name not in row]
        if missing:
            raise AnalysisError(f"Population result row {index} is missing {', '.join(missing)}.")
        source_name = row.get("source_parameter_name", row["parameter_name"])
        if source_name != row["parameter_name"]:
            raise AnalysisError(
                f"Population result row {index} relabels {source_name} as {row['parameter_name']}."
            )
        source_grade = row.get("source_evidence_grade", row["evidence_grade"])
        if source_grade != row["evidence_grade"]:
            raise AnalysisError(
                f"Population result row {index} relabels evidence grade {source_grade} as {row['evidence_grade']}."
            )
        if row["evidence_grade"] not in {"estimated", "predicted"}:
            raise AnalysisError(
                f"Population result row {index} keeps evidence grade estimated or predicted. "
                f"{row['evidence_grade']!r} is rejected."
            )
        payload = {name: row[name] for name in _REQUIRED}
        checked.append(PopulationResultRow.model_validate(payload))
    return PopulationResultTable(
        rows=checked,
        limitations=[
            "The evidence grade is the grade declared by the caller.",
            "Reading this table does not start a subprocess and does not interpret an engine response.",
            "CL and CL/F stay in separate parameter names.",
            "The table is not a dosing recommendation.",
        ],
    )
