"""Noise-free direct-effect Emax, in vitro ADME, and one toxicology endpoint.

Research use only. Hill is fixed at 1. The effect is percent and is paired
with an explicit concentration, not with AUC.
"""

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION, TOX_RESEARCH_LIMITATION
from pharmacoweave.dmpk import AdmeAssayKind, AdmeStatus, InVitroAdmeRecord
from pharmacoweave.pd_models import (
    EffectDimension,
    EffectDirection,
    EffectObservation,
    EffectTimeScale,
    emax_effect,
)
from pharmacoweave.pkpd import ExposureKind, ExposureValue, fit_linked_direct_emax
from pharmacoweave.toxicology import ToxEndpoint, ToxKind, ToxSeverity

E0 = 10.0
EMAX = 80.0
EC50 = 8.0
CONCENTRATIONS = (1.0, 5.0, 10.0, 20.0, 40.0)


def build_effects() -> list[EffectObservation]:
    """Observed effects on the noise-free direct Emax curve."""

    effects: list[EffectObservation] = []
    for index, concentration in enumerate(CONCENTRATIONS, start=1):
        value = emax_effect(concentration=concentration, e0=E0, emax=EMAX, ec50=EC50, hill=1.0)
        effects.append(
            EffectObservation(
                effect_id=f"EFF-{index}",
                endpoint_name="synthetic response",
                matrix="plasma",
                method="synthetic direct effect",
                time_scale=EffectTimeScale.DIRECT,
                dimension=EffectDimension.PERCENT,
                direction=EffectDirection.INCREASE,
                status=ObservationStatus.QUANTIFIED,
                exposure_ref=f"CONC-{index}",
                limitations=[RESEARCH_USE_LIMITATION],
                numeric_value=value,
                unit="percent",
                time_since_exposure=0.0,
                time_since_exposure_unit="h",
            )
        )
    return effects


def build_exposures() -> list[ExposureValue]:
    """Explicit concentrations. These are not AUC values."""

    return [
        ExposureValue(
            source_id=f"CONC-{index}",
            kind=ExposureKind.CONCENTRATION,
            value=concentration,
            unit="mg/L",
        )
        for index, concentration in enumerate(CONCENTRATIONS, start=1)
    ]


def build_adme() -> list[InVitroAdmeRecord]:
    """Three measured in vitro rows and one row that was not measured."""

    shared = {"compound_id": "PW-0001", "limitations": [RESEARCH_USE_LIMITATION]}
    return [
        InVitroAdmeRecord(
            record_id="ADME-PERM",
            assay_kind=AdmeAssayKind.PERMEABILITY,
            system="Caco-2",
            status=AdmeStatus.MEASURED,
            value=12.0e-6,
            unit="cm/s",
            **shared,
        ),
        InVitroAdmeRecord(
            record_id="ADME-MIC",
            assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
            system="human liver microsomes",
            status=AdmeStatus.MEASURED,
            value=25.0,
            unit="uL/min/mg_protein",
            **shared,
        ),
        InVitroAdmeRecord(
            record_id="ADME-HEP",
            assay_kind=AdmeAssayKind.HEPATOCYTE_CLEARANCE,
            system="human hepatocytes",
            status=AdmeStatus.MEASURED,
            value=8.0,
            unit="uL/min/million_cells",
            **shared,
        ),
        InVitroAdmeRecord(
            record_id="ADME-UNSET",
            assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
            system="human liver microsomes",
            status=AdmeStatus.NOT_MEASURED,
            **shared,
        ),
    ]


def build_tox() -> ToxEndpoint:
    """A histopathology finding linked to a dose event. Not a clinical dose."""

    return ToxEndpoint(
        endpoint_id="TOX-1",
        compound_id="PW-0001",
        kind=ToxKind.HISTOPATHOLOGY,
        finding="Minimal hepatocellular hypertrophy in a synthetic example.",
        severity=ToxSeverity.MINIMAL,
        limitations=[TOX_RESEARCH_LIMITATION],
        subject_id="SUBJ-01",
        dose_event_id="DOSE-SUBJ-01",
    )


def main() -> None:
    """Print the recovered direct-effect parameters."""

    _linked, fit = fit_linked_direct_emax(build_effects(), build_exposures(), e0=E0, e0_unit="percent")
    found = {item.name: item.value for item in fit.estimates}
    print(
        f"status={fit.diagnostic.status.value} Emax={found['Emax']:.6g} "
        f"EC50={found['EC50']:.6g} predictions={len(fit.predicted)}"
    )


if __name__ == "__main__":
    main()
