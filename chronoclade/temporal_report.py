"""Temporal-signal assessment, exports and figures."""

from __future__ import annotations

import csv
import math
import re
from html import escape
from io import StringIO
from pathlib import Path


RATE_DETAIL = re.compile(
    r"--rate:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)(?:\s*\+/-\s*([0-9.]+(?:e[+-]?\d+)?))?",
    re.IGNORECASE,
)
R_SQUARED_DETAIL = re.compile(r"--r\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)
CHI_SQUARED_DETAIL = re.compile(r"--chi\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)


def assess_temporal_signal(
    temporal: dict[str, object], *, p_value_threshold: float
) -> dict[str, object]:
    """Return a stable verdict and reader-facing explanation."""

    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    rate = float(observed["rate"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])
    method = str(temporal.get("method", "root_to_tip"))

    if requested == 0:
        return {
            "code": "not_assessed",
            "label": "Temporal signal not assessed",
            "reason": "No date randomisations were requested, so no dated tree was generated.",
            "supported": False,
        }
    if successful != requested:
        return {
            "code": "not_assessed",
            "label": "Temporal signal assessment incomplete",
            "reason": (
                f"Only {successful} of {requested} date randomisations completed successfully. "
                "The temporal-signal gate was not evaluated."
            ),
            "supported": False,
        }
    if rate <= 0:
        return {
            "code": "not_supported",
            "label": "Temporal-signal screen failed",
            "reason": (
                "The estimated evolutionary rate was not positive. Sampling dates do not "
                "support reliable time scaling for this lineage."
            ),
            "supported": False,
        }
    if method == "full_tree":
        if not bool(temporal.get("cr2_passed", False)):
            overlaps = int(temporal.get("cr2_overlapping_randomisations", 0))
            return {
                "code": "not_supported",
                "label": "Full TreeTime randomisation screen failed",
                "reason": (
                    "The observed approximate 95% rate interval overlapped "
                    f"{overlaps} interval(s) from the date-randomised full model fits."
                ),
                "supported": False,
            }
        return {
            "code": "supported",
            "label": "Full TreeTime randomisation screen passed",
            "reason": (
                "The observed approximate 95% rate interval did not overlap any interval "
                "from the date-randomised full TreeTime fits (strict CR2 rule)."
            ),
            "supported": True,
        }

    p_value = float(temporal["p_value_r_squared"])
    if p_value > p_value_threshold:
        return {
            "code": "not_supported",
            "label": "Root-to-tip permutation screen failed",
            "reason": (
                f"The observed root-to-tip fit was not stronger than randomised dates at the "
                f"predefined threshold (p={p_value:.3g}; threshold={p_value_threshold:.3g})."
            ),
            "supported": False,
        }
    return {
        "code": "supported",
        "label": "Root-to-tip permutation screen passed",
        "reason": (
            "The lineage has a positive estimated rate and its observed root-to-tip fit was "
            "stronger than expected after randomly permuting sampling dates."
        ),
        "supported": True,
    }


def observed_metric(temporal: dict[str, object], name: str) -> float:
    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    return float(observed[name])


def format_rate(value: float) -> str:
    return f"{value:.3g} substitutions/site/year"


def _pyplot():
    """Load matplotlib with a non-interactive backend only when figures are written."""

    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = "chronoclade"
    from matplotlib import pyplot

    return pyplot


def _plot_style(axis: object) -> None:
    axis.set_facecolor("#ffffff")
    axis.grid(axis="y", color="#d7d9df", linewidth=0.7, alpha=0.75)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#343842")
    axis.spines["bottom"].set_color("#343842")
    axis.tick_params(colors="#343842", labelsize=9)


def write_randomisation_csv(temporal: dict[str, object], output: Path) -> Path:
    """Write the complete observed and permuted clock results as tidy CSV."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "method",
                "analysis",
                "permutation",
                "rate",
                "rate_std",
                "rate_lower_95",
                "rate_upper_95",
                "r_squared",
                "exceeds_observed",
                "overlaps_observed_rate_interval",
            ],
        )
        writer.writeheader()
        observed_rate = observed_metric(temporal, "rate")
        observed_r_squared = observed_metric(temporal, "r_squared")
        writer.writerow(
            {
                "method": temporal.get("method", "root_to_tip"),
                "analysis": "observed",
                "permutation": "",
                "rate": observed_rate,
                "r_squared": observed_r_squared,
                "exceeds_observed": "",
                "rate_std": temporal["observed"].get("rate_std", ""),
                "rate_lower_95": temporal["observed"].get("rate_lower_95", ""),
                "rate_upper_95": temporal["observed"].get("rate_upper_95", ""),
                "overlaps_observed_rate_interval": "",
            }
        )
        for index, value in enumerate(randomised, start=1):
            r_squared = float(value["r_squared"])
            writer.writerow(
                {
                    "method": temporal.get("method", "root_to_tip"),
                    "analysis": "date_randomisation",
                    "permutation": index,
                    "rate": float(value["rate"]),
                    "r_squared": r_squared,
                    "exceeds_observed": str(r_squared >= observed_r_squared).lower(),
                    "rate_std": value.get("rate_std", ""),
                    "rate_lower_95": value.get("rate_lower_95", ""),
                    "rate_upper_95": value.get("rate_upper_95", ""),
                    "overlaps_observed_rate_interval": (
                        str(
                            max(
                                float(value.get("rate_lower_95", float("inf"))),
                                float(temporal["observed"].get("rate_lower_95", float("-inf"))),
                            )
                            <= min(
                                float(value.get("rate_upper_95", float("-inf"))),
                                float(temporal["observed"].get("rate_upper_95", float("inf"))),
                            )
                        ).lower()
                        if temporal.get("method") == "full_tree"
                        else ""
                    ),
                }
            )
    return output


def write_randomisation_png(temporal: dict[str, object], output: Path) -> Path:
    """Write a publication-ready raster companion to the SVG randomisation plot."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    pyplot = _pyplot()
    figure, axes = pyplot.subplots(1, 2, figsize=(10.6, 4.5), constrained_layout=True)
    panels = (
        (
            axes[0],
            [float(value["r_squared"]) for value in randomised],
            observed_metric(temporal, "r_squared"),
            "Root-to-tip fit",
            "R²",
        ),
        (
            axes[1],
            [float(value["rate"]) for value in randomised],
            observed_metric(temporal, "rate"),
            "Clock rate",
            "Substitutions/site/year",
        ),
    )
    for axis, values, observed, title, label in panels:
        _plot_style(axis)
        if values:
            bins = min(20, max(5, round(math.sqrt(len(values)))))
            axis.hist(values, bins=bins, color="#2855a6", alpha=0.78, edgecolor="white")
        axis.axvline(observed, color="#d63c2f", linewidth=2.4, label="Observed")
        axis.set_title(title, loc="left", fontsize=12, fontweight="bold", color="#17191f")
        axis.set_xlabel(label, color="#343842")
        axis.set_ylabel("Permutations", color="#343842")
        axis.legend(frameon=False, fontsize=9)
    if temporal.get("method") == "full_tree":
        interval_axis = axes[1]
        interval_axis.clear()
        _plot_style(interval_axis)
        for index, value in enumerate(randomised, start=1):
            lower = float(value["rate_lower_95"])
            upper = float(value["rate_upper_95"])
            interval_axis.hlines(index, lower, upper, color="#6285a3", linewidth=1.2, alpha=0.6)
            interval_axis.plot(float(value["rate"]), index, ".", color="#2855a6", markersize=3)
        observed = temporal["observed"]
        observed_y = len(randomised) + 2
        interval_axis.hlines(
            observed_y,
            float(observed["rate_lower_95"]),
            float(observed["rate_upper_95"]),
            color="#d63c2f",
            linewidth=3,
            label="Observed approximate 95% interval",
        )
        interval_axis.plot(float(observed["rate"]), observed_y, "o", color="#d63c2f", markersize=5)
        interval_axis.set_title(
            "Rate-interval separation", loc="left", fontsize=12, fontweight="bold", color="#17191f"
        )
        interval_axis.set_xlabel("Substitutions/site/year", color="#343842")
        interval_axis.set_ylabel("Randomised full fit", color="#343842")
        interval_axis.legend(frameon=False, fontsize=9)
    title = (
        "Full TreeTime tip-date randomisation · strict CR2 "
        + ("passed" if bool(temporal.get("cr2_passed")) else "failed")
        if temporal.get("method") == "full_tree"
        else f"Root-to-tip date randomisation · empirical p={float(temporal['p_value_r_squared']):.3g}"
    )
    figure.suptitle(
        title,
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
        color="#17191f",
    )
    figure.savefig(
        output,
        dpi=220,
        facecolor="#ffffff",
        metadata={"Title": "Date randomisation", "Author": "ChronoClade"},
    )
    pyplot.close(figure)
    return output


def write_root_to_tip_png(source: Path, output: Path) -> Path | None:
    """Render TreeTime's root-to-tip table as a readable PNG scatter plot."""

    if not source.is_file():
        return None
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    points: list[tuple[float, float, str]] = []
    for row in rows:
        normalized = {str(key).strip(): value for key, value in row.items()}
        try:
            points.append(
                (
                    float(str(normalized["date"]).strip()),
                    float(str(normalized["root-to-tip distance"]).strip()),
                    str(normalized.get("name", "")).strip(),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    if not points:
        return None
    dates = [point[0] for point in points]
    distances = [point[1] for point in points]
    mean_x = sum(dates) / len(dates)
    mean_y = sum(distances) / len(distances)
    denominator = sum((value - mean_x) ** 2 for value in dates)
    slope = (
        sum(
            (x_value - mean_x) * (y_value - mean_y)
            for x_value, y_value in zip(dates, distances, strict=True)
        )
        / denominator
        if denominator
        else 0.0
    )
    intercept = mean_y - slope * mean_x
    pyplot = _pyplot()
    figure, axis = pyplot.subplots(figsize=(9.4, 5.4), constrained_layout=True)
    _plot_style(axis)
    axis.scatter(dates, distances, s=34, color="#2855a6", alpha=0.78, linewidth=0)
    lower, upper = min(dates), max(dates)
    axis.plot(
        [lower, upper],
        [intercept + slope * lower, intercept + slope * upper],
        color="#d63c2f",
        linewidth=2.2,
    )
    axis.set_title(
        "Genetic divergence through sampling time",
        loc="left",
        fontsize=15,
        fontweight="bold",
        color="#17191f",
    )
    axis.set_xlabel("Collection date", color="#343842")
    axis.set_ylabel("Root-to-tip distance (substitutions/site)", color="#343842")
    figure.savefig(
        output,
        dpi=220,
        facecolor="#ffffff",
        metadata={"Title": "Root-to-tip regression", "Author": "ChronoClade"},
    )
    pyplot.close(figure)
    return output


def _clock_confidence(directory: Path) -> dict[str, object]:
    """Read TreeTime clock and node-date confidence outputs when available."""

    clock_path = directory / "timetree" / "molecular_clock.txt"
    dates_path = directory / "timetree" / "dates.tsv"
    result: dict[str, object] = {
        "method": "TreeTime maximum likelihood",
        "interval_level": "90% max-posterior region",
    }
    if clock_path.is_file():
        text = clock_path.read_text(encoding="utf-8")
        rate = RATE_DETAIL.search(text)
        r_squared = R_SQUARED_DETAIL.search(text)
        chi_squared = CHI_SQUARED_DETAIL.search(text)
        if rate:
            result["rate"] = float(rate.group(1))
            result["rate_std_dev"] = float(rate.group(2)) if rate.group(2) else None
        if r_squared:
            result["r_squared"] = float(r_squared.group(1))
        if chi_squared:
            result["chi_squared"] = float(chi_squared.group(1))
    node_rows: list[dict[str, object]] = []
    if dates_path.is_file():
        with dates_path.open(encoding="utf-8") as handle:
            reader = csv.reader(
                (line for line in handle if not line.startswith("#")), delimiter="\t"
            )
            for row in reader:
                if len(row) < 5:
                    continue
                try:
                    numeric, lower, upper = map(float, row[2:5])
                except ValueError:
                    continue
                node_rows.append(
                    {
                        "node": row[0],
                        "date": row[1],
                        "numeric_date": numeric,
                        "lower_bound": lower,
                        "upper_bound": upper,
                        "interval_width_years": upper - lower,
                    }
                )
        internal = [row for row in node_rows if str(row["node"]).startswith("NODE_")]
        if internal:
            root = internal[0]
            widths = sorted(float(row["interval_width_years"]) for row in internal)
            result.update(
                {
                    "root_node": root["node"],
                    "root_date": root["date"],
                    "root_numeric_date": root["numeric_date"],
                    "root_lower_bound": root["lower_bound"],
                    "root_upper_bound": root["upper_bound"],
                    "internal_nodes": len(internal),
                    "median_internal_interval_width_years": widths[len(widths) // 2],
                    "maximum_internal_interval_width_years": max(widths),
                }
            )
    result["node_rows"] = node_rows
    return result


def write_timetree_confidence(directory: Path) -> dict[str, object]:
    """Export TreeTime confidence summaries and per-node dates as CSV."""

    confidence = _clock_confidence(directory)
    node_rows = confidence.pop("node_rows", [])
    summary_path = directory / "timetree_confidence.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key, value in confidence.items():
            writer.writerow([key, "" if value is None else value])
    if node_rows:
        node_path = directory / "node_dates.csv"
        with node_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(node_rows[0]))
            writer.writeheader()
            writer.writerows(node_rows)
    return confidence


def write_timetree_figures(directory: Path, confidence: dict[str, object]) -> Path | None:
    """Render a calendar-time phylogeny with visible node-date intervals."""

    source = directory / "timetree" / "timetree.nexus"
    if not source.is_file():
        return None
    from Bio import Phylo

    node_dates = directory / "node_dates.csv"
    if not node_dates.is_file():
        return None
    with node_dates.open(newline="", encoding="utf-8") as handle:
        date_rows = {
            row["node"]: {
                "date": float(row["numeric_date"]),
                "lower": float(row["lower_bound"]),
                "upper": float(row["upper_bound"]),
            }
            for row in csv.DictReader(handle)
        }

    nexus_text = source.read_text(encoding="utf-8")
    tree_line = next(
        line for line in nexus_text.splitlines() if line.lstrip().casefold().startswith("tree ")
    )
    tree = Phylo.read(StringIO(tree_line.split("=", 1)[1].strip()), "newick")
    pyplot = _pyplot()
    tip_count = tree.count_terminals()
    figure, axis = pyplot.subplots(
        figsize=(12.2, max(7.0, min(14.0, tip_count * 0.15))), constrained_layout=True
    )
    terminals = tree.get_terminals()
    y_positions = {clade: float(index) for index, clade in enumerate(reversed(terminals))}

    def y_position(clade: object) -> float:
        if clade in y_positions:
            return y_positions[clade]
        children = list(clade.clades)
        value = sum(y_position(child) for child in children) / len(children)
        y_positions[clade] = value
        return value

    y_position(tree.root)

    def numeric_date(clade: object, parent_date: float | None = None) -> float:
        row = date_rows.get(str(clade.name))
        if row:
            return row["date"]
        if parent_date is None:
            return float(confidence.get("root_numeric_date", 0.0))
        return parent_date + float(clade.branch_length or 0.0)

    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or not clade.name or str(clade.name) not in date_rows:
            continue
        row = date_rows[str(clade.name)]
        axis.hlines(
            y_positions[clade],
            row["lower"],
            row["upper"],
            color="#e2aaa3",
            linewidth=4.2,
            alpha=0.72,
            zorder=1,
        )

    def draw_clade(clade: object, parent_date: float | None = None) -> None:
        x_value = numeric_date(clade, parent_date)
        children = list(clade.clades)
        if children:
            child_y = [y_positions[child] for child in children]
            axis.vlines(
                x_value,
                min(child_y),
                max(child_y),
                color="#17191f",
                linewidth=0.65,
                zorder=2,
            )
        for child in children:
            child_x = numeric_date(child, x_value)
            axis.hlines(
                y_positions[child],
                x_value,
                child_x,
                color="#17191f",
                linewidth=0.7,
                zorder=2,
            )
            draw_clade(child, x_value)

    draw_clade(tree.root)
    tip_x = [numeric_date(clade) for clade in terminals]
    tip_y = [y_positions[clade] for clade in terminals]
    axis.scatter(tip_x, tip_y, s=11, color="#2855a6", zorder=3, label="Sampled genome")
    root_x = numeric_date(tree.root)
    axis.scatter([root_x], [y_positions[tree.root]], s=32, color="#d63c2f", zorder=4, label="Root")
    root_lower = confidence.get("root_lower_bound")
    root_upper = confidence.get("root_upper_bound")
    if root_lower is not None and root_upper is not None:
        axis.annotate(
            f"Root {root_x:.1f}\n90% interval {float(root_lower):.1f}–{float(root_upper):.1f}",
            xy=(root_x, y_positions[tree.root]),
            xytext=(12, 12),
            textcoords="offset points",
            fontsize=9,
            color="#8f241d",
            fontweight="bold",
        )
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="#17191f", linewidth=1, label="Time-scaled branch"),
        Line2D([0], [0], color="#e2aaa3", linewidth=4.2, label="90% node-date interval"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#2855a6",
            markeredgecolor="none",
            label="Sampled genome",
        ),
    ]
    axis.legend(handles=legend_handles, loc="upper left", frameon=False, ncol=3, fontsize=9)
    axis.set_xlabel("Calendar year")
    axis.set_yticks([])
    axis.set_ylabel("")
    axis.set_title(
        "Dated phylogeny with node-time uncertainty",
        loc="left",
        fontsize=15,
        fontweight="bold",
        pad=26,
    )
    axis.grid(axis="x", color="#d7d9df", linewidth=0.7, alpha=0.8)
    axis.margins(x=0.03, y=0.015)
    axis.spines["left"].set_visible(False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_color("#343842")
    output = directory / "timetree_with_confidence.svg"
    figure.savefig(
        output,
        facecolor="white",
        metadata={"Date": None, "Creator": "ChronoClade"},
    )
    figure.savefig(
        directory / "timetree.png",
        dpi=220,
        facecolor="white",
        metadata={
            "Title": "Dated phylogeny with node-time uncertainty",
            "Author": "ChronoClade",
        },
    )
    pyplot.close(figure)
    return output


def _histogram_panel(
    values: list[float],
    *,
    observed: float,
    x: float,
    title: str,
    x_label: str,
    formatter: str,
) -> str:
    panel_width = 430.0
    plot_x = x + 52
    plot_y = 72.0
    plot_width = 340.0
    plot_height = 220.0
    title_text = escape(title)
    label_text = escape(x_label)

    if not values:
        return f"""
        <g>
          <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
          <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
          <text x="{x + panel_width / 2}" y="185" text-anchor="middle" class="empty">No completed permutations</text>
          <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
        </g>
        """

    lower = min([*values, observed])
    upper = max([*values, observed])
    if math.isclose(lower, upper):
        padding = abs(lower) * 0.1 or 1.0
        lower -= padding
        upper += padding
    else:
        padding = (upper - lower) * 0.05
        lower -= padding
        upper += padding

    bin_count = min(20, max(5, round(math.sqrt(len(values)))))
    bin_width = (upper - lower) / bin_count
    counts = [0] * bin_count
    for value in values:
        index = min(bin_count - 1, max(0, int((value - lower) / bin_width)))
        counts[index] += 1
    maximum = max(counts) or 1
    bar_width = plot_width / bin_count
    bars: list[str] = []
    for index, count in enumerate(counts):
        height = plot_height * count / maximum
        bars.append(
            f'<rect x="{plot_x + index * bar_width + 1:.2f}" '
            f'y="{plot_y + plot_height - height:.2f}" width="{max(1, bar_width - 2):.2f}" '
            f'height="{height:.2f}" class="bar"/>'
        )

    observed_x = plot_x + ((observed - lower) / (upper - lower)) * plot_width
    lower_text = format(lower, formatter)
    upper_text = format(upper, formatter)
    observed_text = format(observed, formatter)
    return f"""
    <g>
      <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
      <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
      {"".join(bars)}
      <line x1="{observed_x:.2f}" y1="{plot_y - 8}" x2="{observed_x:.2f}" y2="{plot_y + plot_height}" class="observed"/>
      <text x="{observed_x:.2f}" y="{plot_y - 16}" text-anchor="middle" class="observed-label">Observed {escape(observed_text)}</text>
      <line x1="{plot_x}" y1="{plot_y + plot_height}" x2="{plot_x + plot_width}" y2="{plot_y + plot_height}" class="axis"/>
      <text x="{plot_x}" y="{plot_y + plot_height + 20}" text-anchor="start" class="tick">{escape(lower_text)}</text>
      <text x="{plot_x + plot_width}" y="{plot_y + plot_height + 20}" text-anchor="end" class="tick">{escape(upper_text)}</text>
      <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
    </g>
    """


def write_randomisation_plot(temporal: dict[str, object], output: Path) -> Path:
    """Write an accessible two-panel SVG comparing observed and permuted metrics."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    r_squared = [float(value["r_squared"]) for value in randomised]
    rates = [float(value["rate"]) for value in randomised]
    observed_r_squared = observed_metric(temporal, "r_squared")
    observed_rate = observed_metric(temporal, "rate")
    method = str(temporal.get("method", "root_to_tip"))
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])

    left = _histogram_panel(
        r_squared,
        observed=observed_r_squared,
        x=20,
        title="Root-to-tip fit under randomised dates",
        x_label="R²",
        formatter=".3g",
    )
    if method == "full_tree":
        observed = temporal["observed"]
        all_bounds = [
            bound
            for value in [observed, *randomised]
            for bound in (float(value["rate_lower_95"]), float(value["rate_upper_95"]))
        ]
        lower, upper = min(all_bounds), max(all_bounds)
        span = upper - lower or 1.0

        def x_position(value: float) -> float:
            return 495 + (value - lower) / span * 370

        interval_lines = []
        count = max(1, len(randomised))
        for index, value in enumerate(randomised):
            y = 92 + index / count * 176
            interval_lines.append(
                f'<line x1="{x_position(float(value["rate_lower_95"])):.2f}" y1="{y:.2f}" '
                f'x2="{x_position(float(value["rate_upper_95"])):.2f}" y2="{y:.2f}" '
                'stroke="#6285a3" stroke-width="1.3" opacity=".6"/>'
            )
        observed_y = 286
        right = f"""<g>
          <text x="470" y="28" class="title">Rate-interval separation</text>
          <rect x="485" y="52" width="390" height="252" class="plot-bg"/>
          {''.join(interval_lines)}
          <line x1="{x_position(float(observed['rate_lower_95'])):.2f}" y1="{observed_y}" x2="{x_position(float(observed['rate_upper_95'])):.2f}" y2="{observed_y}" class="observed"/>
          <circle cx="{x_position(observed_rate):.2f}" cy="{observed_y}" r="4" fill="#c6403d"/>
          <text x="485" y="326" class="tick">{lower:.2e}</text><text x="875" y="326" text-anchor="end" class="tick">{upper:.2e}</text>
          <text x="680" y="350" text-anchor="middle" class="axis-label">Substitutions/site/year · observed interval in red</text>
        </g>"""
    else:
        right = _histogram_panel(
            rates,
            observed=observed_rate,
            x=470,
            title="Clock rate under randomised dates",
            x_label="Substitutions/site/year",
            formatter=".2e",
        )
    footer = (
        f"{successful} of {requested} full TreeTime refits completed · strict CR2 "
        + ("passed" if bool(temporal.get("cr2_passed")) else "failed")
        if method == "full_tree"
        else f"{successful} of {requested} permutations completed · root-to-tip randomisation p={float(temporal['p_value_r_squared']):.3g}"
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="920" height="390" viewBox="0 0 920 390" role="img" aria-labelledby="title description">
  <title id="title">Date-randomisation comparison</title>
  <desc id="description">Histograms compare the observed root-to-tip R squared and clock rate with values from permuted sampling dates.</desc>
  <style>
    text {{ font-family: Archivo, sans-serif; fill: #243340; }}
    .title {{ font-size: 16px; font-weight: 700; }}
    .axis-label {{ font-size: 13px; font-weight: 600; }}
    .tick, .empty {{ font-size: 12px; fill: #607080; }}
    .plot-bg {{ fill: #f5f8fa; stroke: #d6e0e6; }}
    .bar {{ fill: #6285a3; opacity: 0.82; }}
    .observed {{ stroke: #c6403d; stroke-width: 3; }}
    .observed-label {{ font-size: 12px; font-weight: 700; fill: #a52f2c; }}
    .axis {{ stroke: #607080; stroke-width: 1; }}
    .footer {{ font-size: 13px; fill: #40515e; }}
  </style>
  {left}
  {right}
  <text x="460" y="372" text-anchor="middle" class="footer">{escape(footer)}</text>
</svg>
"""
    output.write_text(svg, encoding="utf-8")
    return output
