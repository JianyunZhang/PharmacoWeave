"""HTML rendering for the data-validation and research analysis report."""

from datetime import datetime

from jinja2 import Environment, PackageLoader, select_autoescape

from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.core.provenance import OmittedArtifact
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import ValidationReport


def render_validation_report(
    study: PharmacologyStudy,
    report: ValidationReport,
    *,
    run_id: str,
    generated_at: datetime,
    counts: dict[str, int],
    omitted: list[OmittedArtifact],
    version: str,
    schema_version: str,
    nca_parameters: list[NcaParameterSet] | None = None,
    fit_rows: list[dict[str, str]] | None = None,
    emax_rows: list[dict[str, str]] | None = None,
    semilog_svg: str | None = None,
    goodness_of_fit_svg: str | None = None,
    residual_diagnostics_svg: str | None = None,
    weighted_residuals_svg: str | None = None,
    concentration_effect_svg: str | None = None,
) -> str:
    """Render the validation inventory, plus NCA, PK fit, or Emax sections when run."""

    environment = Environment(
        loader=PackageLoader("pharmacoweave.reporting", "templates"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = environment.get_template("validation_report.html.j2")
    description = study.metadata.description or "No study description was provided."
    return template.render(
        study_id=study.metadata.study_id,
        title=study.metadata.title,
        data_class=study.metadata.data_class.value,
        species=study.metadata.species,
        description=description,
        run_id=run_id,
        summary=report.summary(),
        issues=[issue.model_dump() for issue in report.issues],
        counts=counts,
        omitted=omitted,
        version=version,
        schema_version=schema_version,
        generated_at=generated_at.isoformat(),
        nca_rows=[_nca_row(record) for record in nca_parameters or []],
        fit_rows=fit_rows or [],
        emax_rows=emax_rows or [],
        semilog_svg=semilog_svg,
        goodness_of_fit_svg=goodness_of_fit_svg,
        residual_diagnostics_svg=residual_diagnostics_svg,
        weighted_residuals_svg=weighted_residuals_svg,
        concentration_effect_svg=concentration_effect_svg,
    )


def _nca_row(record: NcaParameterSet) -> dict[str, str]:
    return {
        "subject_id": record.subject_id,
        "auc_last": _fmt(record.auc_last),
        "auc_inf": _fmt(record.auc_inf),
        "auc_inf_extrap_fraction": _fmt(record.auc_inf_extrap_fraction),
        "half_life": _fmt(record.half_life),
        "clearance": _fmt(record.clearance),
        "volume": _fmt(record.volume),
        "cmax": _fmt(record.cmax),
        "tmax": _fmt(record.tmax),
    }


def _fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"
