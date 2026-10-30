"""Structured research-use toxicology endpoints.

Importing this module does not produce a clinical dose or a safety claim.
"""

from pharmacoweave.core.roadmap import status
from pharmacoweave.toxicology.records import ToxEndpoint, ToxEndpointSet, ToxKind, ToxSeverity, dosing_recommendation

STATUS = status("toxicology")

__all__ = ["STATUS", "ToxEndpoint", "ToxEndpointSet", "ToxKind", "ToxSeverity", "dosing_recommendation"]
