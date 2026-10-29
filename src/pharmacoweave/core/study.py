"""A pharmacology study and the operations that stay on its boundary."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pydantic import Field, PrivateAttr

if TYPE_CHECKING:
    from pharmacoweave.core.modeling import AnalysisRecords
    from pharmacoweave.core.provenance import ResultBundleManifest
    from pharmacoweave.core.validation import ValidationReport
    from pharmacoweave.pk_models.fit import IvBolusFit

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
from pharmacoweave.core.evidence import EvidenceRecord, IdentityLink, SourceFile


class PharmacologyStudy(Record):
    """One study in memory.

    ``from_files`` reads a directory. ``check`` checks foreign keys, clocks,
    assay limits, and evidence boundaries. ``run_nca`` calculates NCA for the
    locked single-dose IV bolus design. ``export_bundle`` writes a validation
    inventory. ``export_analysis_bundle`` also runs that NCA, fits CL and V,
    and sets ``analysis_executed``.
    """

    metadata: StudyMetadata
    subjects: list[Subject] = Field(default_factory=list)
    compounds: list[Compound] = Field(default_factory=list)
    analytes: list[Analyte] = Field(default_factory=list)
    assays: list[Assay] = Field(default_factory=list)
    dose_events: list[DoseEvent] = Field(default_factory=list)
    biosamples: list[Biosample] = Field(default_factory=list)
    observations: list[Observation] = Field(default_factory=list)
    identity_links: list[IdentityLink] = Field(default_factory=list)
    external_evidence: list[EvidenceRecord] = Field(default_factory=list)
    _source_files: list[SourceFile] = PrivateAttr(default_factory=list)

    @property
    def source_files(self) -> tuple[SourceFile, ...]:
        """Files hashed when this study was read from a directory."""

        return tuple(self._source_files)

    def attach_source_files(self, files: list[SourceFile]) -> None:
        """Record input-file hashes. This does not change scientific fields."""

        self._source_files = list(files)

    @classmethod
    def from_files(
        cls,
        path: str | Path,
        *,
        validate_schema: bool = True,
        mode: str = "research",
    ) -> PharmacologyStudy:
        """Load a study directory.

        With ``validate_schema=True``, relational errors raise
        ``ValidationFailed``. Warnings do not. ``mode='strict'`` promotes
        warnings to errors.
        """

        from pharmacoweave.errors import ValidationFailed
        from pharmacoweave.ingestion.tabular import read_study_directory

        study = read_study_directory(path)
        if validate_schema:
            report = study.check(mode=mode)
            if not report.ok:
                raise ValidationFailed(report)
        return study

    def check(self, *, mode: str = "research") -> ValidationReport:
        """Return a validation report for this study.

        The method is named ``check`` so it does not override Pydantic's
        ``BaseModel.validate`` classmethod.
        """

        from pharmacoweave.core.validation import validate_study

        return validate_study(self, mode=mode)

    def export_bundle(self, path: str | Path, *, run_id: str, mode: str = "research") -> ResultBundleManifest:
        """Write a research-use result bundle. Fails when validation has errors."""

        from pharmacoweave.core.bundle import export_bundle

        return export_bundle(self, path, run_id=run_id, mode=mode)

    def run_nca(self, config: str | Path) -> AnalysisRecords:
        """Calculate NCA for the locked single-dose IV bolus design."""

        from pharmacoweave.pk_nca.runner import run_nca

        return run_nca(self, config)

    def fit_one_compartment_iv_bolus(
        self,
        *,
        error_model: str = "proportional",
        max_nfev: int = 200,
        initial: tuple[float, float] | None = None,
    ) -> IvBolusFit:
        """Estimate CL and V. Non-convergence omits estimates and is not a pass."""

        from pharmacoweave.pk_models.fit import fit_one_compartment_iv_bolus

        return fit_one_compartment_iv_bolus(
            self,
            error_model=error_model,
            max_nfev=max_nfev,
            initial=initial,
        )

    def export_analysis_bundle(
        self,
        path: str | Path,
        *,
        run_id: str,
        nca_config: str | Path,
        mode: str = "research",
        error_model: str = "proportional",
    ) -> ResultBundleManifest:
        """Run the locked IV bolus analysis and write a new result bundle."""

        from pharmacoweave.analysis import export_analysis_bundle

        return export_analysis_bundle(
            self,
            path,
            run_id=run_id,
            nca_config=nca_config,
            mode=mode,
            error_model=error_model,
        )

    def attach_candidates(self, path: str | Path) -> PharmacologyStudy:
        """Attach a TheraWeave evidence file. Predicted ADMET stays predicted."""

        from pharmacoweave.integrations.theraweave import attach_candidate_bundle

        return attach_candidate_bundle(self, path)
