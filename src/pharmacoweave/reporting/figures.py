"""Concentration–time SVG without an extra plotting dependency."""

from pharmacoweave.core.profiles import ConcentrationProfile

_COLORS = ("#0f6f78", "#c46b2c", "#3d4c7a", "#6b4c7a")


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


def _empty(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<text x="24" y="40" font-family="sans-serif" font-size="14">No quantified concentrations.</text>'
        "</svg>"
    )
