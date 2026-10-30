"""SVG figures without an extra plotting dependency.

BLQ concentrations and censored effects are never drawn. Callers pass quantified
points only, and these functions also drop non-positive concentrations instead
of plotting them at zero.
"""

import math
from dataclasses import dataclass
from xml.sax.saxutils import escape

from pharmacoweave.core.profiles import ConcentrationProfile

_COLORS = ("#0f6f78", "#c46b2c", "#3d4c7a", "#6b4c7a")
_BLQ_NOTE = "BLQ concentrations are omitted and are not plotted as zero."
_EFFECT_BLQ_NOTE = "BLQ observations are omitted and are not plotted as zero."


@dataclass(frozen=True)
class PredictedPoint:
    """One quantified observation and its model prediction."""

    series_id: str
    time_h: float
    observed_mg_l: float
    predicted_mg_l: float


@dataclass(frozen=True)
class SemilogOverlay:
    """Optional smooth fitted curve on a semi-log concentration-time figure."""

    series_id: str
    times_h: tuple[float, ...]
    predicted_mg_l: tuple[float, ...]


@dataclass(frozen=True)
class EffectPlotPoint:
    """One quantified concentration-effect pair for PD plotting."""

    series_id: str
    concentration: float
    effect: float


@dataclass(frozen=True)
class EffectCurveOverlay:
    """Fitted direct or sigmoid Emax curve parameters for overlay."""

    series_id: str
    e0: float
    emax: float
    ec50: float
    hill: float = 1.0


def concentration_time_svg(profiles: list[ConcentrationProfile]) -> str:
    """Draw quantified concentrations. Censored points are not plotted as zero."""

    series = [
        (profile.subject_id, [pt for pt in profile.points if pt.concentration_mg_l > 0 and pt.time_h >= 0])
        for profile in profiles
    ]
    series = [(sid, pts) for sid, pts in series if pts]
    width, height = 640, 380
    left, right, top, bottom = 56, 16, 24, 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    if not series:
        return _empty(width, height)
    tmax = max(point.time_h for _sid, pts in series for point in pts)
    cmax = max(point.concentration_mg_l for _sid, pts in series for point in pts)
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
    for index, (subject_id, pts) in enumerate(series):
        color = _COLORS[index % len(_COLORS)]
        coords = " ".join(f"{x_of(point.time_h):.2f},{y_of(point.concentration_mg_l):.2f}" for point in pts)
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{coords}"/>')
        for point in pts:
            cx = f"{x_of(point.time_h):.2f}"
            cy = f"{y_of(point.concentration_mg_l):.2f}"
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="3" fill="{color}"/>')
        legend_y = top + 16 + index * 16
        parts.append(
            f'<text x="{left + 8}" y="{legend_y}" font-family="sans-serif" font-size="12" '
            f'fill="{color}">{escape(subject_id)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def concentration_time_semilog_svg(
    profiles: list[ConcentrationProfile],
    *,
    lloq_mg_l: float | None = None,
    overlays: list[SemilogOverlay] | None = None,
    terminal_point_counts: dict[str, int] | None = None,
) -> str:
    """Draw semi-log (log10) concentration vs time with optional LLOQ, fit overlays, and terminal phase."""

    series = [
        (
            profile.subject_id,
            [
                pt
                for pt in profile.points
                if pt.concentration_mg_l > 0
                and math.isfinite(pt.concentration_mg_l)
                and math.isfinite(pt.time_h)
                and pt.time_h >= 0
            ],
        )
        for profile in profiles
    ]
    series = [(sid, pts) for sid, pts in series if pts]
    width, height = 640, 380
    if not series:
        return _empty(width, height)

    all_c = [pt.concentration_mg_l for _sid, pts in series for pt in pts]
    all_t = [pt.time_h for _sid, pts in series for pt in pts]
    usable_overlays: list[tuple[str, list[tuple[float, float]]]] = []
    for ov in overlays or []:
        pairs = [
            (t, c)
            for t, c in zip(ov.times_h, ov.predicted_mg_l, strict=True)
            if c > 0 and math.isfinite(c) and math.isfinite(t) and t >= 0
        ]
        if pairs:
            usable_overlays.append((ov.series_id, pairs))
            all_c.extend(c for _t, c in pairs)
            all_t.extend(t for t, _c in pairs)

    has_lloq = lloq_mg_l is not None and math.isfinite(lloq_mg_l) and lloq_mg_l > 0
    if has_lloq and lloq_mg_l is not None:
        all_c.append(lloq_mg_l * 0.85)

    min_c = min(all_c)
    max_c = max(all_c)
    d_min = math.floor(math.log10(min_c))
    d_max = math.ceil(math.log10(max_c))
    if d_max <= d_min:
        d_max = d_min + 1
    tmax = max(max(all_t), 1e-6)

    frame = _frame(
        width,
        height,
        x_max=tmax,
        y_min=float(d_min),
        y_max=float(d_max),
        x_label="Time since dose (h)",
        y_label="log10 Concentration (mg/L)",
        title="Semi-log concentration vs time",
    )

    for decade in range(d_min, d_max + 1):
        gy = frame.y_of(float(decade))
        label = f"{10.0**decade:g}"
        frame.parts.append(
            f'<line x1="{frame.left}" y1="{gy:.2f}" x2="{frame.left + frame.plot_w}" y2="{gy:.2f}" '
            'stroke="#e5ecf0" stroke-width="1"/>'
        )
        frame.parts.append(
            f'<text x="{frame.left - 6}" y="{gy + 4:.2f}" text-anchor="end" '
            f'font-family="sans-serif" font-size="11" fill="#486070">{escape(label)}</text>'
        )

    if has_lloq and lloq_mg_l is not None:
        ly = frame.y_of(math.log10(lloq_mg_l))
        frame.parts.append(
            f'<line x1="{frame.left}" y1="{ly:.2f}" x2="{frame.left + frame.plot_w}" y2="{ly:.2f}" '
            'stroke="#b05252" stroke-dasharray="5 3" stroke-width="1.2"/>'
        )
        frame.parts.append(
            f'<text x="{frame.left + frame.plot_w - 4}" y="{ly - 4:.2f}" text-anchor="end" '
            'font-family="sans-serif" font-size="11" fill="#b05252">LLOQ</text>'
        )

    color_by_sid = {sid: _COLORS[idx % len(_COLORS)] for idx, (sid, _pts) in enumerate(series)}
    for ov_sid, pairs in usable_overlays:
        color = color_by_sid.get(ov_sid, _COLORS[0])
        coords = " ".join(f"{frame.x_of(t):.2f},{frame.y_of(math.log10(c)):.2f}" for t, c in pairs)
        frame.parts.append(
            f'<polyline class="fit-overlay" fill="none" stroke="{color}" stroke-width="1.6" '
            f'stroke-dasharray="6 3" points="{coords}"/>'
        )

    term_counts = terminal_point_counts or {}
    for index, (subject_id, pts) in enumerate(series):
        color = _COLORS[index % len(_COLORS)]
        coords = " ".join(
            f"{frame.x_of(pt.time_h):.2f},{frame.y_of(math.log10(pt.concentration_mg_l)):.2f}" for pt in pts
        )
        frame.parts.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="1.8" opacity="0.85" points="{coords}"/>'
        )
        n_term = term_counts.get(subject_id, 0)
        term_set = {pt.observation_id for pt in pts[-n_term:]} if n_term > 0 else set()
        for pt in pts:
            cx = frame.x_of(pt.time_h)
            cy = frame.y_of(math.log10(pt.concentration_mg_l))
            frame.parts.append(f'<circle class="observed-point" cx="{cx:.2f}" cy="{cy:.2f}" r="3.2" fill="{color}"/>')
            if pt.observation_id in term_set:
                frame.parts.append(
                    f'<rect class="terminal-point" x="{cx - 5.5:.2f}" y="{cy - 5.5:.2f}" '
                    f'width="11" height="11" fill="none" stroke="{color}" stroke-width="1.4"/>'
                )
        legend_y = frame.top + 16 + index * 16
        frame.parts.append(
            f'<text x="{frame.left + 8}" y="{legend_y}" font-family="sans-serif" font-size="12" '
            f'fill="{color}">{escape(subject_id)}</text>'
        )

    frame.parts.append("</svg>")
    return "\n".join(frame.parts)


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


def weighted_residual_svg(points: list[PredictedPoint]) -> str:
    """Relative weighted residual ((observed - predicted) / predicted) versus time."""

    usable = _quantified(points)
    width, height = 640, 380
    if not usable:
        return _empty(width, height)
    wres = [(point.observed_mg_l - point.predicted_mg_l) / point.predicted_mg_l for point in usable]
    tmax = max(max(point.time_h for point in usable), 1e-6)
    span = max(max(wres) - min(wres), 1e-6)
    pad = 0.1 * span
    y_min = min(min(wres) - pad, -0.05)
    y_max = max(max(wres) + pad, 0.05)
    frame = _frame(
        width,
        height,
        x_max=tmax,
        y_min=y_min,
        y_max=y_max,
        x_label="Time since dose (h)",
        y_label="Relative residual ((obs - pred) / pred)",
        title="Weighted residual vs time",
    )
    zero_y = frame.y_of(0.0)
    frame.parts.append(
        f'<line x1="{frame.left}" y1="{zero_y:.2f}" x2="{frame.left + frame.plot_w}" y2="{zero_y:.2f}" '
        'stroke="#9aa7b2" stroke-dasharray="4 3"/>'
    )
    _scatter(
        frame,
        usable,
        lambda point: (point.time_h, (point.observed_mg_l - point.predicted_mg_l) / point.predicted_mg_l),
    )
    frame.parts.append("</svg>")
    return "\n".join(frame.parts)


def concentration_effect_svg(
    points: list[EffectPlotPoint],
    *,
    overlays: list[EffectCurveOverlay] | None = None,
    concentration_unit: str = "mg/L",
    effect_unit: str = "effect",
) -> str:
    """Draw quantified concentration-effect pairs and optional fitted direct/sigmoid Emax curves."""

    usable = [
        pt
        for pt in points
        if pt.concentration > 0 and math.isfinite(pt.concentration) and math.isfinite(pt.effect)
    ]
    width, height = 640, 380
    if not usable:
        return _empty_effect(width, height)

    valid_overlays = [
        ov
        for ov in (overlays or [])
        if ov.ec50 > 0
        and ov.hill > 0
        and math.isfinite(ov.e0)
        and math.isfinite(ov.emax)
        and math.isfinite(ov.ec50)
        and math.isfinite(ov.hill)
    ]
    x_max = max(pt.concentration for pt in usable) * 1.05
    x_max = max(x_max, 1e-6)
    y_candidates = [pt.effect for pt in usable]
    for ov in valid_overlays:
        y_candidates.extend([ov.e0, ov.e0 + ov.emax])
    y_lo, y_hi = min(y_candidates), max(y_candidates)
    span = max(y_hi - y_lo, 1e-6)
    pad = 0.08 * span
    frame = _frame(
        width,
        height,
        x_max=x_max,
        y_min=y_lo - pad,
        y_max=y_hi + pad,
        x_label=f"Concentration ({concentration_unit})",
        y_label=f"Effect ({effect_unit})",
        title="Concentration vs effect",
        desc=_EFFECT_BLQ_NOTE,
    )

    series_ids = list(dict.fromkeys(pt.series_id for pt in usable))
    color_by_sid = {sid: _COLORS[idx % len(_COLORS)] for idx, sid in enumerate(series_ids)}

    for ov in valid_overlays:
        color = color_by_sid.get(ov.series_id, _COLORS[0])
        base_y = frame.y_of(ov.e0)
        frame.parts.append(
            f'<line x1="{frame.left}" y1="{base_y:.2f}" x2="{frame.left + frame.plot_w}" y2="{base_y:.2f}" '
            'stroke="#9aa7b2" stroke-dasharray="3 3"/>'
        )
        if 0 < ov.ec50 <= x_max:
            ec50_x = frame.x_of(ov.ec50)
            frame.parts.append(
                f'<line x1="{ec50_x:.2f}" y1="{frame.top}" x2="{ec50_x:.2f}" y2="{frame.top + frame.plot_h}" '
                f'stroke="{color}" stroke-dasharray="4 3" opacity="0.65"/>'
            )
            frame.parts.append(
                f'<text x="{ec50_x + 4:.2f}" y="{frame.top + frame.plot_h - 6:.2f}" '
                f'font-family="sans-serif" font-size="11" fill="{color}">EC50</text>'
            )
        curve_pts: list[str] = []
        n_grid = 64
        for step in range(n_grid + 1):
            c_val = x_max * step / n_grid
            eff_val = _emax_eval(c_val, ov.e0, ov.emax, ov.ec50, ov.hill)
            curve_pts.append(f"{frame.x_of(c_val):.2f},{frame.y_of(eff_val):.2f}")
        frame.parts.append(
            f'<polyline class="emax-curve" fill="none" stroke="{color}" stroke-width="2" '
            f'points="{" ".join(curve_pts)}"/>'
        )

    for index, sid in enumerate(series_ids):
        color = _COLORS[index % len(_COLORS)]
        selected = [pt for pt in usable if pt.series_id == sid]
        for pt in selected:
            cx = f"{frame.x_of(pt.concentration):.2f}"
            cy = f"{frame.y_of(pt.effect):.2f}"
            frame.parts.append(f'<circle cx="{cx}" cy="{cy}" r="3.5" fill="{color}"/>')
        legend_y = frame.top + 16 + index * 16
        frame.parts.append(
            f'<text x="{frame.left + 8}" y="{legend_y}" font-family="sans-serif" font-size="12" '
            f'fill="{color}">{escape(sid)}</text>'
        )

    frame.parts.append("</svg>")
    return "\n".join(frame.parts)


def _emax_eval(c_val: float, e0: float, emax: float, ec50: float, hill: float) -> float:
    if c_val <= 0.0:
        return e0
    z = hill * (math.log(c_val) - math.log(ec50))
    if z >= 40.0:
        g = 1.0
    elif z <= -40.0:
        g = 0.0
    elif z >= 0.0:
        g = 1.0 / (1.0 + math.exp(-z))
    else:
        ez = math.exp(z)
        g = ez / (1.0 + ez)
    return e0 + emax * g


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
    desc: str = _BLQ_NOTE,
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
            f"<desc>{escape(desc)} Research use only.</desc>",
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


def _empty_effect(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        "<title>No quantified effects</title>"
        f"<desc>{escape(_EFFECT_BLQ_NOTE)}</desc>"
        '<text x="24" y="40" font-family="sans-serif" font-size="14">No quantified concentration-effect pairs.</text>'
        "</svg>"
    )
