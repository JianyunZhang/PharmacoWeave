"""CSV, TSV, and Parquet IO for a study directory.

The column tuples in ``pharmacoweave.core.tables`` are the file contract.
Absolute times stay ISO-8601 strings. A blank cell is an absent value, not zero.
"""

import csv
import hashlib
from pathlib import Path
from typing import Any, cast

import pyarrow as pa
import pyarrow.parquet as pq

from pharmacoweave.core.entities import (
    Analyte,
    Assay,
    Biosample,
    Compound,
    DoseEvent,
    Observation,
    Record,
    StudyMetadata,
    Subject,
)
from pharmacoweave.core.evidence import IdentityLink, SourceFile
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.tables import (
    BOOL_COLUMNS,
    FLOAT_COLUMNS,
    OPTIONAL_TABLES,
    REQUIRED_TABLES,
    TABLE_COLUMNS,
    arrow_schema,
)
from pharmacoweave.errors import IngestError

_MODELS: dict[str, type[Record]] = {
    "subjects": Subject,
    "compounds": Compound,
    "analytes": Analyte,
    "assays": Assay,
    "dosing_events": DoseEvent,
    "biosamples": Biosample,
    "observations": Observation,
    "identity_map": IdentityLink,
}
_SUFFIXES = (".csv", ".tsv", ".parquet")


def read_study_directory(path: str | Path) -> PharmacologyStudy:
    """Read ``study_metadata.json`` and the study tables from ``path``."""

    directory = Path(path)
    if not directory.is_dir():
        raise IngestError(f"{directory} is not a study directory.")
    metadata_path = directory / "study_metadata.json"
    if not metadata_path.is_file():
        raise IngestError(f"{directory} is missing study_metadata.json.")
    try:
        metadata = StudyMetadata.model_validate_json(metadata_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise IngestError(f"study_metadata.json does not match StudyMetadata: {exc}") from exc

    located: dict[str, Path] = {"study_metadata": metadata_path}
    for stem in REQUIRED_TABLES:
        found = _locate(directory, stem, required=True)
        assert found is not None
        located[stem] = found
    identity_path = _locate(directory, "identity_map", required=False)
    if identity_path is not None:
        located["identity_map"] = identity_path

    loaded = {stem: _read_table(file_path, stem) for stem, file_path in located.items() if stem != "study_metadata"}
    study = PharmacologyStudy(
        metadata=metadata,
        subjects=cast(list[Subject], loaded["subjects"]),
        compounds=cast(list[Compound], loaded["compounds"]),
        analytes=cast(list[Analyte], loaded["analytes"]),
        assays=cast(list[Assay], loaded["assays"]),
        dose_events=cast(list[DoseEvent], loaded["dosing_events"]),
        biosamples=cast(list[Biosample], loaded["biosamples"]),
        observations=cast(list[Observation], loaded["observations"]),
        identity_links=cast(list[IdentityLink], loaded.get("identity_map", [])),
    )
    study.attach_source_files(
        [
            SourceFile(relative_path=file_path.name, sha256=_sha256(file_path), role=role)
            for role, file_path in located.items()
        ]
    )
    return study


def write_study_directory(study: PharmacologyStudy, path: str | Path, *, file_format: str = "csv") -> None:
    """Write a study directory. ``file_format`` is ``csv``, ``tsv``, or ``parquet``."""

    if file_format not in {"csv", "tsv", "parquet"}:
        raise IngestError("file_format must be csv, tsv, or parquet.")
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    metadata_path = directory / "study_metadata.json"
    metadata_path.write_text(study.metadata.model_dump_json(indent=2) + "\n", encoding="utf-8")
    rows = {
        "subjects": [subject_row(item) for item in study.subjects],
        "compounds": [compound_row(item) for item in study.compounds],
        "analytes": [analyte_row(item) for item in study.analytes],
        "assays": [assay_row(item) for item in study.assays],
        "dosing_events": [dose_row(item) for item in study.dose_events],
        "biosamples": [sample_row(item) for item in study.biosamples],
        "observations": [observation_row(item) for item in study.observations],
        "identity_map": [identity_row(item) for item in study.identity_links],
    }
    for stem, records in rows.items():
        _write_table(directory, stem, TABLE_COLUMNS[stem], records, file_format)


def write_named_table(path: str | Path, stem: str, rows: list[dict[str, object]]) -> None:
    """Write one Parquet table using the column contract for ``stem``."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    _write_parquet(destination, TABLE_COLUMNS[stem], rows)


def observations_frame(study: PharmacologyStudy) -> object:
    """Return observations as a flat pandas DataFrame."""

    import pandas as pd

    return pd.DataFrame(
        [observation_row(item) for item in study.observations],
        columns=list(TABLE_COLUMNS["observations"]),
    )


def subject_row(subject: Subject) -> dict[str, object]:
    weight = subject.body_weight
    return {
        "subject_id": subject.subject_id,
        "study_id": subject.study_id,
        "group_id": subject.group_id,
        "species": subject.species,
        "sex": subject.sex.value,
        "body_weight_value": None if weight is None else weight.value,
        "body_weight_unit": None if weight is None else weight.unit,
    }


def compound_row(compound: Compound) -> dict[str, object]:
    return compound.model_dump(mode="json")


def analyte_row(analyte: Analyte) -> dict[str, object]:
    return analyte.model_dump(mode="json")


def assay_row(assay: Assay) -> dict[str, object]:
    return assay.model_dump(mode="json")


def dose_row(dose: DoseEvent) -> dict[str, object]:
    return dose.model_dump(mode="json")


def sample_row(sample: Biosample) -> dict[str, object]:
    return sample.model_dump(mode="json")


def observation_row(observation: Observation) -> dict[str, object]:
    payload = observation.model_dump(mode="json")
    payload["qc_flags"] = "|".join(observation.qc_flags)
    return payload


def identity_row(link: IdentityLink) -> dict[str, object]:
    return link.model_dump(mode="json")


def _locate(directory: Path, stem: str, *, required: bool) -> Path | None:
    found = [directory / f"{stem}{suffix}" for suffix in _SUFFIXES if (directory / f"{stem}{suffix}").is_file()]
    if len(found) > 1:
        names = ", ".join(path.name for path in found)
        raise IngestError(f"Multiple files supply {stem}: {names}.")
    if not found:
        if required or stem in REQUIRED_TABLES:
            raise IngestError(f"Missing {stem}.csv, {stem}.tsv, or {stem}.parquet.")
        if stem not in OPTIONAL_TABLES and stem not in REQUIRED_TABLES:
            raise IngestError(f"Unknown table {stem}.")
        return None
    return found[0]


def _read_table(path: Path, stem: str) -> list[Any]:
    columns = TABLE_COLUMNS[stem]
    model = _MODELS[stem]
    if path.suffix == ".parquet":
        table = pq.read_table(path)
        if tuple(table.column_names) != columns:
            raise IngestError(f"{path.name} columns {tuple(table.column_names)} do not match the contract.")
        raw_rows = table.to_pylist()
    else:
        delimiter = "\t" if path.suffix == ".tsv" else ","
        raw_rows = _read_delimited(path, columns, delimiter)
    built: list[Any] = []
    for lineno, raw in enumerate(raw_rows, start=2):
        payload = _coerce_row(raw, columns, path.name, lineno)
        if stem == "subjects":
            value = payload.pop("body_weight_value")
            unit = payload.pop("body_weight_unit")
            payload["body_weight"] = None if value is None and unit is None else {"value": value, "unit": unit}
        try:
            built.append(model.model_validate(payload))
        except Exception as exc:
            raise IngestError(f"{path.name} row {lineno}: {exc}") from exc
    return built


def _read_delimited(path: Path, columns: tuple[str, ...], delimiter: str) -> list[dict[str, str | None]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise IngestError(f"{path.name} has no header.")
        names = tuple(name.strip() for name in reader.fieldnames)
        if names != columns:
            raise IngestError(f"{path.name} columns {names} do not match {columns}.")
        return [{key: row[key] for key in columns} for row in reader]


def _coerce_row(
    row: dict[str, object],
    columns: tuple[str, ...],
    filename: str,
    lineno: int,
) -> dict[str, object]:
    prepared: dict[str, object] = {}
    for column in columns:
        value = _clean(row.get(column))
        try:
            if column in FLOAT_COLUMNS and isinstance(value, str):
                value = float(value)
            elif column in BOOL_COLUMNS and isinstance(value, str):
                value = _parse_bool(value)
            elif column == "qc_flags":
                if value is None:
                    value = []
                elif isinstance(value, str):
                    value = [part for part in value.split("|") if part]
        except ValueError as exc:
            raise IngestError(f"{filename} row {lineno} column {column}: {exc}") from exc
        if isinstance(value, float) and not _finite(value):
            raise IngestError(f"{filename} row {lineno} column {column} is not finite.")
        prepared[column] = value
    return prepared


def _write_table(
    directory: Path,
    stem: str,
    columns: tuple[str, ...],
    rows: list[dict[str, object]],
    file_format: str,
) -> None:
    if file_format == "parquet":
        _write_parquet(directory / f"{stem}.parquet", columns, rows)
        return
    delimiter = "\t" if file_format == "tsv" else ","
    _write_delimited(directory / f"{stem}.{file_format}", columns, rows, delimiter)


def _write_delimited(
    path: Path,
    columns: tuple[str, ...],
    rows: list[dict[str, object]],
    delimiter: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _cell(row.get(column)) for column in columns})


def _write_parquet(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    normalized: list[dict[str, object]] = []
    for row in rows:
        item: dict[str, object] = {}
        for column in columns:
            value = row.get(column)
            if isinstance(value, list):
                value = "|".join(str(part) for part in value)
            item[column] = value
        normalized.append(item)
    schema = arrow_schema(columns)
    table = pa.Table.from_pylist(normalized, schema=schema)
    pq.write_table(table, path)


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return format(value, ".10g")
    if isinstance(value, list):
        return "|".join(str(part) for part in value)
    return str(value)


def _clean(value: object) -> object:
    if isinstance(value, str):
        stripped = value.strip()
        return None if stripped == "" else stripped
    return value


def _parse_bool(value: str) -> bool:
    text = value.strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"Expected true or false, found {value!r}.")


def _finite(value: float) -> bool:
    return value == value and value not in {float("inf"), float("-inf")}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
