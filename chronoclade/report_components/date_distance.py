"""Shared figures for descriptive distance-versus-date diagnostics."""

from __future__ import annotations

import math
from pathlib import Path
import textwrap


def _finite(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _location(point: dict) -> str:
    return str(point.get("country") or point.get("location") or "Unknown")


def _label(point: dict) -> str:
    return str(point.get("label") or point.get("sample_id") or "Unnamed sample")


def _style(axis) -> None:
    axis.set_facecolor("#ffffff")
    axis.grid(axis="y", color="#d7d9df", linewidth=0.7, alpha=0.75)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=9)


def _point(axis, point: dict, x: float, y: float, colours: dict[str, str]) -> None:
    colour = colours[_location(point)]
    is_input = point.get("role") == "input"
    axis.plot(
        x,
        y,
        marker="o" if is_input else "^",
        linestyle="none",
        markersize=6,
        markerfacecolor=colour,
        markeredgecolor="#17191f" if is_input else colour,
        markeredgewidth=1.1 if is_input else 0.6,
        zorder=3,
    )


def write_date_distance_figure(
    diagnostic: dict,
    samples: list[dict],
    output: Path,
    *,
    title: str,
    distance_label: str,
    colours: dict[str, str],
) -> tuple[Path, Path]:
    """Render supplied points and fits without fitting or inferring any science.

    Date intervals represent collection-date precision, never confidence limits.
    Samples with finite distances but no accepted date retain a separate panel.
    Matplotlib is loaded only when an artifact is requested.
    """
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import ScalarFormatter

    points = diagnostic.get("points", [])
    distance_field = diagnostic.get("distance_field", "distance")
    dated = [
        point
        for point in points
        if _finite(point.get("midpoint_year")) is not None
        and _finite(point.get(distance_field)) is not None
    ]
    dated_ids = {point.get("sample_id") for point in dated}
    undated = [
        sample
        for sample in samples
        if sample.get("sample_id") not in dated_ids
        and _finite(sample.get(distance_field)) is not None
    ]
    rows = 3 if undated else 2
    ratios = [2.7, 1.1, 1.0] if undated else [2.7, 1.1]
    locations = sorted({_location(point) for point in dated + undated})
    legend_rows = math.ceil((len(locations) + 3) / 3)
    legend_height = max(0, legend_rows - 2) * 0.32
    with matplotlib.rc_context({"svg.hashsalt": "chronoclade", "svg.fonttype": "none"}):
        fig, axes = plt.subplots(
            rows,
            1,
            figsize=(10.2, (7.7 if undated else 6.7) + legend_height),
            gridspec_kw={"height_ratios": ratios},
            layout="constrained",
        )
        main, residual = axes[:2]
        for axis in axes:
            _style(axis)
            formatter = ScalarFormatter(useOffset=False)
            axis.xaxis.set_major_formatter(formatter)
        for point in dated:
            x = float(point["midpoint_year"])
            y = float(point[distance_field])
            low = _finite(point.get("year_min"))
            high = _finite(point.get("year_max"))
            colour = colours[_location(point)]
            if low is not None and high is not None:
                main.hlines(y, low, high, color=colour, linewidth=1.2, alpha=0.65)
            _point(main, point, x, y, colours)
            if len(dated) <= 12:
                label = _label(point)
                years = [float(item["midpoint_year"]) for item in dated]
                near_right = x >= min(years) + 0.75 * (max(years) - min(years))
                main.annotate(
                    label if len(label) <= 32 else label[:29] + "…",
                    (x, y),
                    xytext=(-5 if near_right else 5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    ha="right" if near_right else "left",
                )
            value = _finite(point.get("residual"))
            if value is not None:
                if low is not None and high is not None:
                    residual.hlines(value, low, high, color=colour, linewidth=1.2, alpha=0.65)
                _point(residual, point, x, value, colours)
        fitted = diagnostic.get("status") == "fitted"
        if fitted and dated:
            ends = sorted(dated, key=lambda point: float(point["midpoint_year"]))
            # Predictions were computed by the scientific diagnostic, not by the renderer.
            line_points = [point for point in ends if _finite(point.get("predicted")) is not None]
            if line_points:
                main.plot(
                    [float(point["midpoint_year"]) for point in line_points],
                    [float(point["predicted"]) for point in line_points],
                    color="#17191f",
                    linewidth=1.6,
                    zorder=2,
                )
            residual.axhline(0, color="#17191f", linewidth=1)
        else:
            reason = str(diagnostic.get("reason") or "A regression could not be fitted.")
            residual.text(
                0.02,
                0.5,
                textwrap.fill(reason, 95),
                transform=residual.transAxes,
                ha="left",
                va="center",
                fontsize=10,
                color="#5f6470",
            )
            residual.set_yticks([])
        if not dated:
            main.text(
                0.5,
                0.5,
                "No eligible dated points",
                transform=main.transAxes,
                ha="center",
                va="center",
                fontsize=12,
                color="#5f6470",
            )
            main.set_xticks([])
            residual.set_xticks([])
        else:
            residual.set_xlim(main.get_xlim())
        main.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=13)
        main.set_ylabel(distance_label, fontsize=10)
        main.set_xlabel("Collection date (decimal year; interval midpoint)", fontsize=10)
        residual.set_ylabel("Residual", fontsize=10)
        residual.set_title("Residuals (" + distance_label + ")", loc="left", fontsize=9)
        residual.set_xlabel("Collection date (decimal year; interval midpoint)", fontsize=10)
        if undated:
            missing = axes[2]
            for index, sample in enumerate(undated, start=1):
                _point(missing, sample, index, float(sample[distance_field]), colours)
            missing.set_ylabel("Distance", fontsize=10)
            missing.set_title("Undated samples (" + distance_label + ")", loc="left", fontsize=9)
            missing.set_xlabel("Samples without an eligible date (no date assigned)", fontsize=10)
            if len(undated) <= 12:
                missing.set_xticks(range(1, len(undated) + 1))
                missing.set_xticklabels(
                    [
                        _label(point) if len(_label(point)) <= 24 else _label(point)[:21] + "…"
                        for point in undated
                    ],
                    fontsize=8,
                )
            else:
                missing.set_xticks([])
            missing.set_xlim(0.5, len(undated) + 0.5)
        legend = [
            Line2D(
                [],
                [],
                color=colours[location],
                marker="s",
                linestyle="none",
                markersize=7,
                label=textwrap.fill(location, 26),
            )
            for location in locations
        ]
        legend.extend(
            [
                Line2D(
                    [],
                    [],
                    marker="o",
                    markerfacecolor="#ffffff",
                    markeredgecolor="#17191f",
                    linestyle="none",
                    label="Input · circle, dark outline",
                ),
                Line2D(
                    [],
                    [],
                    marker="^",
                    color="#5f6470",
                    linestyle="none",
                    label="Context · triangle",
                ),
            ]
        )
        if fitted:
            legend.append(Line2D([], [], color="#17191f", label="Descriptive fitted line"))
        # Put all entries below the panels so numerous locations cannot cover observations.
        fig.legend(
            handles=legend,
            loc="outside lower center",
            ncol=min(3, len(legend)),
            frameon=False,
            fontsize=9,
        )
        svg = output.with_suffix(".svg")
        png = output.with_suffix(".png")
        fig.savefig(svg, facecolor="white", metadata={"Title": title, "Creator": "ChronoClade"})
        fig.savefig(
            png, dpi=200, facecolor="white", metadata={"Title": title, "Author": "ChronoClade"}
        )
        plt.close(fig)
    return svg, png
