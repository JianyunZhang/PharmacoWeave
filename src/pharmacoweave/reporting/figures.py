"""SVG figures without an extra plotting dependency.

BLQ concentrations are never drawn. Callers pass quantified points only, and
these functions also drop a non-positive observation instead of plotting it
at zero.
"""

import math
from dataclasses import dataclass
from xml.sax.saxutils import escape

from pharmacoweave.core.profiles import ConcentrationProfile

_COLORS = ("#0f6f78", "#c46b2c", "#3d4c7a", "#6b4c7a")
_BLQ_NOTE = "BLQ concentrations are omitted and are not plotted as zero."


def concentration_time_svg(profiles: list[ConcentrationProfile]) -> str:
    """Draw quantified concentrations. Censored points are not plotted as zero."""

    series = [profile for profile in profiles if profile.points]
    width, height = 640, 380
    left, right, top, bottom = 56, 16, 24, 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    if not series:
        return _empty(width, height)
    tmax = max(point.time_h for profile in series for point in profile.points)
    cmax = max(point.concentration_mg_l for profile in series for point in profile.points)
    tmax = max(tmax, 1e-6)
    cmax = max(cmax, 1e-6)

    def x_of(time_h: float) -> float:
        return left + plot_w * time_h / tmax

    def y_of(concentration: float) -> float:
        return top + plot_h * (1.0 - concentration / cmax)

    x_label = f'<text x="{left + plot_w / 2}" y="{height - 12}" text-anchor="middle" '
    x_label += 'font-family="sans-serif" font-size="13">Time since dose (h)</text>'
    y_label = f'<text x="16" y="{top + plot_h / 2}" transform="rotate(-90 16 {top + plot_h / 2})" '
    y_label += 'text-anchor="middle" font-family="sans-serif" font-size="13">Concentration (mg/L)</text>'
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#243140"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#243140"/>',
        x_label,
        y_label,
    ]
    for index, profile in enumerate(series):
        color = _COLORS[index % len(_COLORS)]
        coords = " ".join(f"{x_of(point.time_h):.2f},{y_of(point.concentration_mg_l):.2f}" for point in profile.points)
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{coords}"/>')
        for point in profile.points:
            cx = f"{x_of(point.time_h):.2f}"
            cy = f"{y_of(point.concentration_mg_l):.2f}"
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="3" fill="{color}"/>')
        legend_y = top + 16 + index * 16
        parts.append(
            f'<text x="{left + 8}" y="{legend_y}" font-family="sans-serif" font-size="12" '
            f'fill="{color}">{profile.subject_id}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


@dataclass(frozen=True)
class PredictedPoint:
    """One quantified observation and its model prediction."""

    series_id: str
    time_h: float
    observed_mg_l: float
    predicted_mg_l: float


def goodness_of_fit_svg(points: list[PredictedPoint]) -> str:
    """Observed versus predicted, with the identity line. Non-positive observations are omitted."""

    usable = _quantified(points)
    width, height = 640, 380
    if not usable:
        return _empty(width, height)
    limit = max(max(point.observed_mg_l, point.predicted_mg_l) for point in usable)
    limit = max(limit * 1.05, 1e-6)
    frame = _frame(
        width,
        height,
        x_max=limit,
        y_min=0.0,
        y_max=limit,
        x_label="Observed (mg/L)",
        y_label="Predicted (mg/L)",
        title="Observed vs predicted",
    )
    x0, y0 = frame.x_of(0.0), frame.y_of(0.0)
    x1, y1 = frame.x_of(limit), frame.y_of(limit)
    frame.parts.append(
        f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
        'stroke="#9aa7b2" stroke-dasharray="4 3"/>'
    )
    _scatter(frame, usable, lambda point: (point.observed_mg_l, point.predicted_mg_l))
    frame.parts.append("</svg>")
    return "\n".join(frame.parts)


def residual_diagnostics_svg(points: list[PredictedPoint]) -> str:
    """Residual (observed minus predicted) versus time. Non-positive observations are omitted."""

    usable = _quantified(points)
    width, height = 640, 380
    if not usable:
        return _empty(width, height)
    residuals = [point.observed_mg_l - point.predicted_mg_l for point in usable]
    tmax = max(point.time_h for point in usable)
    tmax = max(tmax, 1e-6)
    span = max(max(residuals) - min(residuals), 1e-6)
    pad = 0.1 * span
    y_min = min(min(residuals) - pad, 0.0)
    y_max = max(max(residuals) + pad, 0.0)
    frame = _frame(
        width,
        height,
        x_max=tmax,
        y_min=y_min,
        y_max=y_max,
        x_label="Time since dose (h)",
        y_label="Residual (mg/L)",
        title="Residual vs time",
    )
    zero_y = frame.y_of(0.0)
    frame.parts.append(
        f'<line x1="{frame.left}" y1="{zero_y:.2f}" x2="{frame.left + frame.plot_w}" y2="{zero_y:.2f}" '
        'stroke="#9aa7b2" stroke-dasharray="4 3"/>'
    )
    _scatter(frame, usable, lambda point: (point.time_h, point.observed_mg_l - point.predicted_mg_l))
    frame.parts.append("</svg>")
    return "\n".join(frame.parts)


@dataclass
class _Frame:
    parts: list[str]
    left: int
    top: int
    plot_w: int
    plot_h: int
    x_max: float
    y_min: float
    y_max: float

    def x_of(self, value: float) -> float:
        return self.left + self.plot_w * value / self.x_max

    def y_of(self, value: float) -> float:
        span = self.y_max - self.y_min
        return self.top + self.plot_h * (self.y_max - value) / span


def _frame(
    width: int,
    height: int,
    *,
    x_max: float,
    y_min: float,
    y_max: float,
    x_label: str,
    y_label: str,
    title: str,
) -> _Frame:
    left, right, top, bottom = 72, 16, 36, 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    frame = _Frame([], left, top, plot_w, plot_h, x_max, y_min, y_max)
    x_text = f'<text x="{left + plot_w / 2}" y="{height - 12}" text-anchor="middle" '
    x_text += f'font-family="sans-serif" font-size="13">{escape(x_label)}</text>'
    y_text = f'<text x="18" y="{top + plot_h / 2}" transform="rotate(-90 18 {top + plot_h / 2})" '
    y_text += f'text-anchor="middle" font-family="sans-serif" font-size="13">{escape(y_label)}</text>'
    frame.parts.extend(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">',
            f"<title>{escape(title)}</title>",
            f"<desc>{escape(_BLQ_NOTE)} Research use only.</desc>",
            '<rect width="100%" height="100%" fill="#ffffff"/>',
            f'<text x="{left}" y="22" font-family="sans-serif" font-size="14">{escape(title)}</text>',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#243140"/>',
            f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#243140"/>',
            x_text,
            y_text,
        ]
    )
    return frame


def _scatter(frame: _Frame, points: list[PredictedPoint], coordinates) -> None:
    series_ids = list(dict.fromkeys(point.series_id for point in points))
    for index, series_id in enumerate(series_ids):
        color = _COLORS[index % len(_COLORS)]
        selected = [point for point in points if point.series_id == series_id]
        for point in selected:
            x_value, y_value = coordinates(point)
            cx = f"{frame.x_of(x_value):.2f}"
            cy = f"{frame.y_of(y_value):.2f}"
            frame.parts.append(f'<circle cx="{cx}" cy="{cy}" r="3.5" fill="{color}"/>')
        legend_y = frame.top + 16 + index * 16
        frame.parts.append(
            f'<text x="{frame.left + 8}" y="{legend_y}" font-family="sans-serif" font-size="12" '
            f'fill="{color}">{escape(series_id)}</text>'
        )


def _quantified(points: list[PredictedPoint]) -> list[PredictedPoint]:
    return [
        point
        for point in points
        if point.observed_mg_l > 0
        and point.predicted_mg_l > 0
        and math.isfinite(point.observed_mg_l)
        and math.isfinite(point.predicted_mg_l)
        and math.isfinite(point.time_h)
        and point.time_h >= 0
    ]


def _empty(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        "<title>No quantified concentrations</title>"
        f"<desc>{escape(_BLQ_NOTE)}</desc>"
        '<text x="24" y="40" font-family="sans-serif" font-size="14">No quantified concentrations.</text>'
        "</svg>"
    )
