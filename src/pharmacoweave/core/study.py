"""A pharmacology study and the operations that stay on its boundary."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import Field, PrivateAttr

if TYPE_CHECKING:
    from pharmacoweave.core.modeling import NcaParameterSet
    from pharmacoweave.core.provenance import ResultBundleManifest
    from pharmacoweave.core.validation import ValidationReport

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
    UrineCollection,
)
from pharmacoweave.core.evidence import EvidenceRecord, IdentityLink, SourceFile


class PharmacologyStudy(Record):
    """One study in memory.

    ``from_files`` reads a directory. ``check`` checks foreign keys, clocks,
    assay limits, and evidence boundaries. ``run_nca`` runs single-dose IV
    bolus NCA.
    """

    metadata: StudyMetadata
    subjects: list[Subject] = Field(default_factory=list)
    compounds: list[Compound] = Field(default_factory=list)
    analytes: list[Analyte] = Field(default_factory=list)
    assays: list[Assay] = Field(default_factory=list)
    dose_events: list[DoseEvent] = Field(default_factory=list)
    biosamples: list[Biosample] = Field(default_factory=list)
    urine_collections: list[UrineCollection] = Field(default_factory=list)
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

    @classmethod
    def from_mapping(
        cls,
        mapping_path: str | Path,
        *,
        base_dir: str | Path | None = None,
        validate_schema: bool = True,
        mode: str = "research",
    ) -> PharmacologyStudy:
        """Load a study from raw laboratory files using an ingest_mapping.yaml specification."""

        from pharmacoweave.ingestion.mapping import ingest_study_from_mapping

        return ingest_study_from_mapping(
            mapping_path,
            base_dir=base_dir,
            mode=mode,
            validate_schema=validate_schema,
        )

    def check(self, *, mode: str = "research") -> ValidationReport:
        """Return a validation report for this study.

        The method is named ``check`` so it does not override Pydantic's
        ``BaseModel.validate`` classmethod.
        """

        from pharmacoweave.core.validation import validate_study

        return validate_study(self, mode=mode)

    def export_bundle(
        self,
        path: str | Path,
        *,
        run_id: str,
        mode: str = "research",
        nca_parameters: list[NcaParameterSet] | None = None,
        fits: Sequence[Any] | None = None,
        emax_fits: Sequence[Any] | None = None,
        adme_records: Sequence[Any] | None = None,
        tox_endpoints: Sequence[Any] | None = None,
        blq_policy: Any | None = None,
        ivive_predictions: Sequence[Any] | None = None,
        prediction_audit: Any | None = None,
        exposure_simulation: Any | None = None,
        urine_results: Sequence[Any] | None = None,
    ) -> ResultBundleManifest:
        """Write a research-use result bundle. Fails when validation has errors."""

        from pharmacoweave.core.bundle import export_bundle

        return export_bundle(
            self,
            path,
            run_id=run_id,
            mode=mode,
            nca_parameters=nca_parameters,
            fits=fits,
            emax_fits=emax_fits,
            adme_records=adme_records,
            tox_endpoints=tox_endpoints,
            blq_policy=blq_policy,
            ivive_predictions=ivive_predictions,
            prediction_audit=prediction_audit,
            exposure_simulation=exposure_simulation,
            urine_results=urine_results,
        )

    def run_nca(self, config: str | Path | Any) -> list[NcaParameterSet]:
        """Run single-dose rich-serial NCA. The study must already validate."""

        from pharmacoweave.pk_nca.runner import run_nca

        return run_nca(self, config)

    def attach_candidates(self, path: str | Path) -> PharmacologyStudy:
        """Attach a TheraWeave evidence file. Predicted ADMET stays predicted."""

        from pharmacoweave.integrations.theraweave import attach_candidate_bundle

        return attach_candidate_bundle(self, path)

    def bridge_effects(
        self, *, endpoint_overrides: dict[str, dict[str, Any]] | None = None
    ) -> list[Any]:
        """Convert study effect observations into typed EffectObservation records."""

        from pharmacoweave.core.observation_bridge import bridge_study_effect_observations

        return bridge_study_effect_observations(self, endpoint_overrides=endpoint_overrides)

    def bridge_toxicology(
        self, *, endpoint_overrides: dict[str, dict[str, Any]] | None = None
    ) -> list[Any]:
        """Convert study toxicity observations into typed ToxEndpoint records."""

        from pharmacoweave.core.observation_bridge import bridge_study_toxicology_endpoints

        return bridge_study_toxicology_endpoints(self, endpoint_overrides=endpoint_overrides)
