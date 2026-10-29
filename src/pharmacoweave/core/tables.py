"""Flat column contracts shared by CSV, TSV, and Parquet.

Absolute times are ISO-8601 strings in every table so the timezone offset
survives language boundaries. ``qc_flags`` is a ``|``-separated string.
"""

import pyarrow as pa

SUBJECT_COLUMNS = (
    "subject_id",
    "study_id",
    "group_id",
    "species",
    "sex",
    "body_weight_value",
    "body_weight_unit",
)
COMPOUND_COLUMNS = (
    "compound_id",
    "name",
    "molecular_weight",
    "molecular_weight_unit",
    "structure_ref",
    "formulation_ref",
    "theraweave_candidate_id",
)
ANALYTE_COLUMNS = (
    "analyte_id",
    "compound_id",
    "name",
    "analyte_role",
    "matrix",
)
ASSAY_COLUMNS = (
    "assay_id",
    "method",
    "analyte_id",
    "matrix",
    "calibration_id",
    "quantitation",
    "lloq",
    "uloq",
    "limit_unit",
    "qc_status",
)
DOSE_COLUMNS = (
    "dose_event_id",
    "subject_id",
    "compound_id",
    "absolute_time",
    "dose_value",
    "dose_unit",
    "route",
    "formulation",
    "infusion_duration",
    "infusion_duration_unit",
    "body_weight_normalized_value",
    "body_weight_normalized_unit",
)
SAMPLE_COLUMNS = (
    "sample_id",
    "subject_id",
    "collection_time_absolute",
    "time_since_dose",
    "time_since_dose_unit",
    "dose_event_id",
    "matrix",
    "collection_site",
    "timebase_id",
    "is_predose",
)
OBSERVATION_COLUMNS = (
    "observation_id",
    "sample_id",
    "analyte_id",
    "observation_type",
    "status",
    "raw_value",
    "reported_value",
    "numeric_value",
    "unit",
    "assay_id",
    "qc_flags",
    "exclusion_reason",
    "replicate_id",
    "assay_batch_id",
    "analyte_form",
)
IDENTITY_COLUMNS = (
    "link_id",
    "source_project",
    "source_study_id",
    "source_subject_id",
    "canonical_subject_id",
    "match_method",
    "review_status",
    "notes",
)

FLOAT_COLUMNS = frozenset(
    {
        "body_weight_value",
        "molecular_weight",
        "lloq",
        "uloq",
        "dose_value",
        "infusion_duration",
        "body_weight_normalized_value",
        "time_since_dose",
        "numeric_value",
    }
)
BOOL_COLUMNS = frozenset({"is_predose"})

REQUIRED_TABLES = (
    "subjects",
    "compounds",
    "analytes",
    "assays",
    "dosing_events",
    "biosamples",
    "observations",
)
OPTIONAL_TABLES = ("identity_map",)
TABLE_COLUMNS = {
    "subjects": SUBJECT_COLUMNS,
    "compounds": COMPOUND_COLUMNS,
    "analytes": ANALYTE_COLUMNS,
    "assays": ASSAY_COLUMNS,
    "dosing_events": DOSE_COLUMNS,
    "biosamples": SAMPLE_COLUMNS,
    "observations": OBSERVATION_COLUMNS,
    "identity_map": IDENTITY_COLUMNS,
}


def arrow_schema(columns: tuple[str, ...]) -> pa.Schema:
    """Build a nullable Arrow schema for one table contract."""

    fields: list[pa.Field] = []
    for name in columns:
        if name in FLOAT_COLUMNS:
            datatype: pa.DataType = pa.float64()
        elif name in BOOL_COLUMNS:
            datatype = pa.bool_()
        else:
            datatype = pa.string()
        fields.append(pa.field(name, datatype, nullable=name not in BOOL_COLUMNS))
    return pa.schema(fields)
