"""Shared offline ordination views; dataset markers stay distinct in every view."""
from __future__ import annotations

from collections import Counter
from html import escape
from pathlib import Path
import re

import numpy as np

from chronoclade.metadata_dates import sample_date_interval


_INPUT_ORIGINS = {"input", "local", "query", "focal"}


def write_ordination_views(output, prefix, coordinates, records, *,
                           title="Distance ordination", axis_labels=("Axis 1", "Axis 2"),
                           groups=None):
    """Write consistent dataset/country/date/group views using identical coordinates.

    Collection dates are midpoint decimal years of their recorded intervals.
    Missing dates remain visible as grey points, never zero-valued dates.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from chronoclade.report_components.styles import report_styles

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", prefix):
        raise ValueError("Ordination prefix must be a safe artifact name")
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (len(records), 2) or not np.isfinite(coordinates).all():
        raise ValueError("Ordination requires two finite coordinates per record")
    ids = [row["sample_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Ordination sample IDs must be unique")
    roles = ["input" if row.get("role", row.get("origin")) in _INPUT_ORIGINS else "context"
             for row in records]
    markers = {"input": "^", "context": "o"}
    names = {"input": "Input genomes", "context": "Public comparisons"}
    dates = [sample_date_interval(row) for row in records]
    years = np.array([value["midpoint_year"]
                      if value["status"] == "valid" else np.nan for value in dates])
    categories = {
        "dataset": [names[role] for role in roles],
        "country": [str(row.get("country") or "Unknown country") for row in records],
        "host": [str(row.get("host") or "Unknown host") for row in records],
    }
    if groups is not None:
        categories["group"] = [str(groups.get(ident, "Ungrouped")) for ident in ids]
    paths, views = {}, []
    for view in [*categories, "date"]:
        fig, ax = plt.subplots(figsize=(9, 6))
        if view == "date":
            finite = np.isfinite(years)
            lower, upper = (float(years[finite].min()), float(years[finite].max())) if finite.any() else (0, 1)
            norm = plt.Normalize(lower, upper if upper > lower else lower + 1)
            for role, marker in markers.items():
                selected = np.array([r == role for r in roles])
                known, missing = selected & finite, selected & ~finite
                if known.any():
                    scatter = ax.scatter(*coordinates[known].T, c=years[known], cmap="viridis",
                                         norm=norm, marker=marker, alpha=.8, s=35,
                                         edgecolors="#222", linewidths=.3)
                if missing.any():
                    ax.scatter(*coordinates[missing].T, c="#999", marker=marker, s=35,
                               edgecolors="#222", linewidths=.3)
            if finite.any():
                fig.colorbar(scatter, ax=ax, label="Collection year (interval midpoint)")
            handles = [Line2D([], [], color="#999", marker="o", linestyle="None", label="Date unknown")]
        else:
            values = categories[view]
            unique = sorted(set(values))
            palette = plt.get_cmap("tab20", max(1, len(unique)))
            colours = {value: palette(i) for i, value in enumerate(unique)}
            counts = Counter(values)
            shown = sorted(unique, key=lambda value: (-counts[value], value))[:20]
            handles = [Line2D([], [], color=colours[value], marker="o", linestyle="None", label=value)
                       for value in shown]
            # One scatter per role avoids quadratic scans and thousands of artists.
            for role, marker in markers.items():
                selected = [i for i, observed in enumerate(roles) if observed == role]
                if selected:
                    ax.scatter(*coordinates[selected].T,
                               c=[colours[values[i]] for i in selected], marker=marker,
                               alpha=.8, s=35, edgecolors="#222", linewidths=.3)
            if len(unique) > len(shown):
                ax.text(0, -.15,
                        f"Legend shows {len(shown)} of {len(unique)} categories by frequency; all genomes are plotted. "
                        "Full labels are in the downloaded tables.",
                        transform=ax.transAxes, fontsize=8, wrap=True)
        handles += [Line2D([], [], color="#222", marker=marker, linestyle="None", label=names[role])
                    for role, marker in markers.items() if role in roles]
        if view == "date":
            ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01),
                      ncol=3, fontsize=8, frameon=False)
        else:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                      fontsize=8, frameon=False)
        ax.set(xlabel=axis_labels[0], ylabel=axis_labels[1], title=title)
        if view == "date":
            ax.set_title(title, pad=30)
        fig.tight_layout()
        filename = f"{prefix}_pcoa" + ("" if view == "dataset" else f"_{view}") + ".svg"
        fig.savefig(output / filename, bbox_inches="tight")
        plt.close(fig)
        key = "pcoa_figure" if view == "dataset" else f"pcoa_{view}_figure"
        paths[key] = str(output / filename)
        views.append((view, filename))
    radios = "".join(f'<label><input type="radio" name="view" value="{view}"'
                     f'{" checked" if i == 0 else ""}>{escape(view.title())}</label> '
                     for i, (view, _) in enumerate(views))
    figures = "".join(f'<figure data-view="{view}"{" hidden" if i else ""}>'
                      f'<img src="{escape(filename, quote=True)}" alt="{escape(title, quote=True)}: {view}"'
                      ' style="width:100%;height:auto"></figure>'
                      for i, (view, filename) in enumerate(views))
    document = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                f'<title>{escape(title)}</title><style>{report_styles()} '
                'body{background:white;padding:8px}figure{margin:0}label{margin-right:18px} '
                'fieldset{border:0;border-bottom:1px solid #ccc}</style></head><body>'
                f'<fieldset><legend>Colour by · marker shape shows input or comparison</legend>{radios}</fieldset>'
                f'{figures}<p class="muted">All views use the same coordinates. Missing metadata remains visible.'
                ' Dates use collection-interval midpoints.</p>'
                '<script>document.querySelectorAll("input[name=view]").forEach(r=>r.addEventListener("change",()=>{'
                'document.querySelectorAll("figure[data-view]").forEach(f=>f.hidden=f.dataset.view!==r.value);'
                '}));</script></body></html>')
    viewer = output / f"{prefix}_ordination.html"
    viewer.write_text(document)
    paths["pcoa_views_html"] = str(viewer)
    return paths


def ordination_viewer(cohorts, directory):
    """Embed only existing contained report assets, with an offline fallback link."""
    directory = Path(directory).resolve()
    frames = []
    for cohort in cohorts or []:
        reference = cohort.get("pcoa_views_html")
        if not reference:
            continue
        path = Path(reference)
        path = (directory / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.is_file() or not path.is_relative_to(directory):
            continue
        relative = escape(path.relative_to(directory).as_posix(), quote=True)
        frames.append(f'<iframe src="{relative}" title="Ordination metadata views" '
                      'style="width:100%;height:760px;border:0" loading="lazy"></iframe>'
                      f'<p><a href="{relative}">Open ordination views</a></p>')
    return "".join(frames)
