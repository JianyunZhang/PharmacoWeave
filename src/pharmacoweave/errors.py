"""Errors raised by PharmacoWeave.

``UnitError`` is also a ``ValueError`` so record constructors surface unit
problems through Pydantic. The other errors are operational and stay out of
that path.
"""


class PharmacoWeaveError(Exception):
    """Base type for errors raised by this package."""


class UnitError(PharmacoWeaveError, ValueError):
    """A unit is unknown or is not convertible to the requested unit."""


class IngestError(PharmacoWeaveError):
    """A study directory or table could not be read into the data model."""


class ValidationFailed(PharmacoWeaveError):
    """A study or bundle failed relational validation."""

    def __init__(self, report: object) -> None:
        self.report = report
        super().__init__(str(report))


class ContractViolation(PharmacoWeaveError):
    """An evidence payload crosses a project boundary it is not allowed to cross."""


class AnalysisError(PharmacoWeaveError):
    """An analysis could not be completed for the locked design."""


class PlannedCapabilityError(PharmacoWeaveError):
    """The data contract for this capability exists, and the calculator does not."""

    def __init__(self, module: str, phase: str, detail: str) -> None:
        self.module = module
        self.phase = phase
        super().__init__(f"{module} ({phase}) is not executable in this version. {detail}")
