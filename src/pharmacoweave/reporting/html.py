"""HTML rendering for the data-validation report."""

from datetime import datetime

from jinja2 import Environment, PackageLoader, select_autoescape

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
    analysis_executed: bool = False,
) -> str:
    """Render the validation inventory. This is not a PK analysis report."""

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
        analysis_executed=analysis_executed,
    )
