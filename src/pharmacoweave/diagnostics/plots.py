"""SVG figures for the IV bolus analysis path.

Figures are drawn from numbers the caller already calculated. A BLQ sample is
marked at the LLOQ when that limit is known, and the marker does not carry a
concentration. Missing LLOQ values are labeled on the time axis instead of
being drawn at zero.
"""

import html
import math
from dataclasses import dataclass

_WIDTH = 760
_HEIGHT = 460
_LEFT = 72
_RIGHT = 24
_TOP = 72
_BOTTOM = 56
_COLORS = ("#0f6f78", "#b45309", "#1d4e89", "#6b3fa0", "#0f766e", "#9f1239")


@dataclass(frozen=True)
class PlotPoint:
    """One sample to draw. ``concentration`` is absent for a censored point."""

    time_h: float
    status: str
    concentration: float | None = None


@dataclass(frozen=True)
class PlotSeries:
    """One subject profile, plus an optional fitted curve in the same unit."""

    name: str
    unit: str
    lloq: float | None
    points: tuple[PlotPoint, ...]
    prediction: tuple[tuple[float, float], ...] = ()


def concentration_time_svg(series: tuple[PlotSeries, ...] | list[PlotSeries], *, title: str) -> str:
    """Draw observed concentrations, censored BLQ markers, and optional fits."""

    quantified = [
        point.concentration
        for item in series
        for point in item.points
        if point.concentration is not None and point.concentration > 0 and point.status == "quantified"
    ]
    predicted = [value for item in series for _time, value in item.prediction if value > 0]
    lloqs = [item.lloq for item in series if item.lloq is not None and item.lloq > 0]
    positive = [*quantified, *predicted, *lloqs]
    times = [point.time_h for item in series for point in item.points]
    times.extend(time for item in series for time, _value in item.prediction)
    if not positive or not times:
        return _document(title, "concentration-time", "<text x='72' y='120'>No positive concentrations to draw.</text>")
    y_min = min(positive)
    y_max = max(positive)
    if y_min == y_max:
        y_min *= 0.8
        y_max *= 1.2
    else:
        log_min = math.log(y_min)
        log_max = math.log(y_max)
        pad = 0.08 * (log_max - log_min)
        y_min = math.exp(log_min - pad)
        y_max = math.exp(log_max + pad)
    x_min = min(0.0, min(times))
    x_max = max(times)
    if x_min == x_max:
        x_max = x_min + 1.0
    else:
        x_max += 0.04 * (x_max - x_min)

    def sx(time_h: float) -> float:
        return _LEFT + (time_h - x_min) / (x_max - x_min) * (_WIDTH - _LEFT - _RIGHT)

    def sy(value: float) -> float:
        return _TOP + (math.log(y_max) - math.log(value)) / (math.log(y_max) - math.log(y_min)) * (
            _HEIGHT - _TOP - _BOTTOM
        )

    parts = [
        _axis_frame(),
        f"<text x='{(_LEFT + _WIDTH - _RIGHT) / 2}' y='{_HEIGHT - 16}' text-anchor='middle' "
        f"font-size='12'>{_esc('Time (h)')}</text>",
        f"<text x='18' y='{(_TOP + _HEIGHT - _BOTTOM) / 2}' text-anchor='middle' font-size='12' "
        f"transform='rotate(-90 18 {(_TOP + _HEIGHT - _BOTTOM) / 2})'>"
        f"{_esc(_concentration_label(series))}</text>",
        "<text x='72' y='62' font-size='12' fill='#486070'>Log concentration scale. "
        "BLQ markers are not measured concentrations.</text>",
    ]
    for index, item in enumerate(series):
        color = _COLORS[index % len(_COLORS)]
        if item.prediction:
            polyline = " ".join(f"{sx(time):.2f},{sy(value):.2f}" for time, value in item.prediction if value > 0)
            parts.append(
                f"<polyline data-series='prediction' data-name='{_esc(item.name)}' fill='none' "
                f"stroke='{color}' stroke-width='1.5' points='{polyline}'/>"
            )
        if item.lloq is not None and item.lloq > 0:
            y_lloq = sy(item.lloq)
            parts.append(
                f"<line data-role='lloq' x1='{_LEFT}' x2='{_WIDTH - _RIGHT}' y1='{y_lloq:.2f}' "
                f"y2='{y_lloq:.2f}' stroke='{color}' stroke-dasharray='4 4' stroke-width='1'/>"
            )
        for point in item.points:
            if point.status == "quantified" and point.concentration is not None and point.concentration > 0:
                parts.append(
                    f"<circle data-status='quantified' data-time='{_num(point.time_h)}' "
                    f"data-conc='{_num(point.concentration)}' cx='{sx(point.time_h):.2f}' "
                    f"cy='{sy(point.concentration):.2f}' r='3.5' fill='{color}'/>"
                )
                continue
            if point.status != "BLQ":
                continue
            if item.lloq is not None and item.lloq > 0:
                x_pos = sx(point.time_h)
                y_pos = sy(item.lloq)
                triangle = (
                    f"{x_pos:.2f},{y_pos - 5:.2f} {x_pos - 4.5:.2f},{y_pos + 4:.2f} "
                    f"{x_pos + 4.5:.2f},{y_pos + 4:.2f}"
                )
                parts.append(
                    f"<polygon data-status='BLQ' data-time='{_num(point.time_h)}' data-conc='' "
                    f"points='{triangle}' fill='none' stroke='{color}' stroke-width='1.4'/>"
                )
            else:
                parts.append(
                    f"<text data-status='BLQ' data-time='{_num(point.time_h)}' data-conc='' "
                    f"x='{sx(point.time_h):.2f}' y='{_HEIGHT - _BOTTOM + 14}' font-size='10'>{_esc('BLQ')}</text>"
                )
        parts.append(
            f"<text x='{_WIDTH - 180}' y='{88 + 16 * index}' font-size='12' fill='{color}'>"
            f"{_esc(item.name)}</text>"
        )
    return _document(title, "concentration-time", "\n".join(parts))


def goodness_of_fit_svg(pairs: list[tuple[float, float]], *, unit: str, title: str) -> str:
    """Draw observed versus predicted concentrations and the identity line."""

    values = [value for pair in pairs for value in pair if value > 0]
    if not values:
        return _document(title, "goodness-of-fit", "<text x='72' y='120'>No paired predictions to draw.</text>")
    lower = min(values) * 0.9
    upper = max(values) * 1.1
    if lower == upper:
        lower *= 0.8
        upper *= 1.2

    def sx(value: float) -> float:
        return _LEFT + (value - lower) / (upper - lower) * (_WIDTH - _LEFT - _RIGHT)

    def sy(value: float) -> float:
        return _TOP + (upper - value) / (upper - lower) * (_HEIGHT - _TOP - _BOTTOM)

    line = f"{sx(lower):.2f},{sy(lower):.2f} {sx(upper):.2f},{sy(upper):.2f}"
    circles = "\n".join(
        f"<circle data-status='paired' data-observed='{_num(observed)}' data-predicted='{_num(predicted)}' "
        f"cx='{sx(predicted):.2f}' cy='{sy(observed):.2f}' r='3.5' fill='#0f6f78'/>"
        for observed, predicted in pairs
    )
    body = "\n".join(
        [
            _axis_frame(),
            f"<polyline fill='none' stroke='#94a3b8' stroke-width='1' points='{line}'/>",
            circles,
            f"<text x='{(_LEFT + _WIDTH - _RIGHT) / 2}' y='{_HEIGHT - 16}' text-anchor='middle' font-size='12'>"
            f"{_esc(f'Predicted ({unit})')}</text>",
            f"<text x='18' y='{(_TOP + _HEIGHT - _BOTTOM) / 2}' text-anchor='middle' font-size='12' "
            f"transform='rotate(-90 18 {(_TOP + _HEIGHT - _BOTTOM) / 2})'>{_esc(f'Observed ({unit})')}</text>",
        ]
    )
    return _document(title, "goodness-of-fit", body)


def residual_svg(points: list[tuple[float, float]], *, title: str) -> str:
    """Draw weighted residuals against time, with a zero reference line."""

    if not points:
        return _document(title, "residual-diagnostics", "<text x='72' y='120'>No residuals to draw.</text>")
    times = [time for time, _residual in points]
    residuals = [residual for _time, residual in points]
    x_min = min(0.0, min(times))
    x_max = max(times) if max(times) > x_min else x_min + 1.0
    span = max(abs(min(residuals)), abs(max(residuals)), 1e-6)
    y_min, y_max = -span * 1.2, span * 1.2

    def sx(time_h: float) -> float:
        return _LEFT + (time_h - x_min) / (x_max - x_min) * (_WIDTH - _LEFT - _RIGHT)

    def sy(value: float) -> float:
        return _TOP + (y_max - value) / (y_max - y_min) * (_HEIGHT - _TOP - _BOTTOM)

    circles = "\n".join(
        f"<circle data-status='residual' data-time='{_num(time_h)}' data-residual='{_num(residual)}' "
        f"cx='{sx(time_h):.2f}' cy='{sy(residual):.2f}' r='3.5' fill='#b45309'/>"
        for time_h, residual in points
    )
    body = "\n".join(
        [
            _axis_frame(),
            f"<line x1='{_LEFT}' x2='{_WIDTH - _RIGHT}' y1='{sy(0):.2f}' y2='{sy(0):.2f}' "
            "stroke='#94a3b8' stroke-width='1'/>",
            circles,
            f"<text x='{(_LEFT + _WIDTH - _RIGHT) / 2}' y='{_HEIGHT - 16}' text-anchor='middle' font-size='12'>"
            "Time (h)</text>",
            f"<text x='18' y='{(_TOP + _HEIGHT - _BOTTOM) / 2}' text-anchor='middle' font-size='12' "
            f"transform='rotate(-90 18 {(_TOP + _HEIGHT - _BOTTOM) / 2})'>Weighted residual</text>",
        ]
    )
    return _document(title, "residual-diagnostics", body)


def _concentration_label(series: list[PlotSeries] | tuple[PlotSeries, ...]) -> str:
    units = sorted({item.unit for item in series if item.unit})
    if len(units) == 1:
        return f"Concentration ({units[0]})"
    return "Concentration"


def _document(title: str, role: str, body: str) -> str:
    return (
        "<?xml version='1.0' encoding='UTF-8'?>\n"
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{_WIDTH}' height='{_HEIGHT}' "
        f"viewBox='0 0 {_WIDTH} {_HEIGHT}' data-role='{role}'>\n"
        f"<rect width='100%' height='100%' fill='#ffffff'/>\n"
        f"<text x='24' y='28' font-family='sans-serif' font-size='16'>{_esc(title)}</text>\n"
        "<text x='24' y='46' font-family='sans-serif' font-size='12' fill='#0f6f78'>"
        "Research use only. Not a dosing recommendation.</text>\n"
        f"<g font-family='sans-serif' fill='#1c2834'>\n{body}\n</g>\n"
        "</svg>\n"
    )


def _axis_frame() -> str:
    return (
        f"<line x1='{_LEFT}' y1='{_TOP}' x2='{_LEFT}' y2='{_HEIGHT - _BOTTOM}' stroke='#1c2834'/>"
        f"<line x1='{_LEFT}' y1='{_HEIGHT - _BOTTOM}' x2='{_WIDTH - _RIGHT}' y2='{_HEIGHT - _BOTTOM}' "
        "stroke='#1c2834'/>"
    )


def _esc(value: str) -> str:
    return html.escape(value, quote=True)


def _num(value: float) -> str:
    return format(value, ".10g")
