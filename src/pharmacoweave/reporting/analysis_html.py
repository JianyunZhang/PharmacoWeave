"""HTML for a research-use PK analysis report."""

from datetime import datetime

from jinja2 import Environment, PackageLoader, select_autoescape

from pharmacoweave.core.modeling import AnalysisRecords
from pharmacoweave.core.provenance import OmittedArtifact


def render_analysis_report(
    records: AnalysisRecords,
    *,
    study_id: str,
    title: str,
    run_id: str,
    generated_at: datetime,
    omitted: list[OmittedArtifact],
    version: str,
    schema_version: str,
) -> str:
    """Render NCA rows, fit status, and omissions. This is not a regulatory report."""

    environment = Environment(
        loader=PackageLoader("pharmacoweave.reporting", "templates"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = environment.get_template("analysis_report.html.j2")
    return template.render(
        study_id=study_id,
        title=title,
        run_id=run_id,
        nca_rows=[item.model_dump(mode="json") for item in records.nca_parameters],
        fits=[item.model_dump(mode="json") for item in records.model_fits],
        estimates=[item.model_dump(mode="json") for item in records.parameter_estimates],
        omitted=omitted,
        version=version,
        schema_version=schema_version,
        generated_at=generated_at.isoformat(),
    )
