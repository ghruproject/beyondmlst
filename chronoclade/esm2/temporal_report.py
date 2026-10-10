"""Offline ESM2 and conventional cgMLST date-diagnostic report."""

from __future__ import annotations

import csv
import json
from html import escape
from pathlib import Path

from chronoclade.location_network.colours import country_palette
from chronoclade.report_components.date_distance import write_date_distance_figure
from chronoclade.report_components.styles import report_styles


def _text(value: object) -> str:
    if value is None or value == "":
        return "Not reported"
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _number(value: object, *, missing: str = "Not fitted") -> str:
    return missing if value is None else f"{float(value):.4g}"


def _records(value: object) -> list[dict]:
    return [item for item in value if isinstance(item, dict)] if isinstance(value, list) else []


def _mapping(value: object) -> dict:
    return value if isinstance(value, dict) else {}


def _sample_rows(source: dict, diagnostic: dict) -> list[dict]:
    """Keep excluded/undated samples and attach accepted-date fitted values."""
    points = _records(diagnostic.get("points"))
    samples = _records(source.get("samples"))
    if not samples:
        return points
    accepted = {point.get("sample_id"): point for point in points}
    rows = [{**sample, **accepted.get(sample.get("sample_id"), {})} for sample in samples]
    ids = {row.get("sample_id") for row in rows}
    return rows + [point for point in points if point.get("sample_id") not in ids]


def _write_csv(path: Path, rows: list[dict]) -> None:
    preferred = [
        "panel",
        "sample_id",
        "label",
        "role",
        "country",
        "location",
        "date_lower",
        "date_upper",
        "precision",
        "midpoint_year",
        "distance",
        "predicted",
        "residual",
        "status",
        "reason",
    ]
    keys = {key for row in rows for key in row}
    fields = [key for key in preferred if key in keys]
    fields += sorted(keys.difference(fields))
    # Even an empty diagnostic has a usable, named table rather than a zero-byte download.
    if not fields:
        fields = preferred
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )


def _table(rows: list[tuple[str, object]], *, label: str) -> str:
    cells = "".join(
        f'<tr><th scope="row">{escape(name)}</th><td>{escape(_text(value))}</td></tr>'
        for name, value in rows
    )
    return (
        f'<div class="table-scroll compact-table" role="region" '
        f'aria-label="{escape(label, quote=True)}" tabindex="0">'
        f"<table><tbody>{cells}</tbody></table></div>"
    )


def _download(name: str, description: str) -> str:
    return (
        f'<li class="download"><a href="{escape(name, quote=True)}" download>'
        f"<span>{escape(name)}</span><small>{escape(description)}</small>"
        f"<b>{escape(Path(name).suffix[1:].upper())}</b></a></li>"
    )


def _samples_table(rows: list[dict], distance_field: str, units: str) -> str:
    if not rows:
        return "<p>No sample observations were available.</p>"
    body = []
    for row in rows:
        interval_data = _mapping(row.get("date_interval"))
        date_low = row.get("start") or row.get("date_start") or interval_data.get("start")
        date_high = row.get("end") or row.get("date_end") or interval_data.get("end")
        interval = _text(date_low)
        if date_high and date_high != date_low:
            interval += " to " + str(date_high)
        label = row.get("label") or row.get("sample_id")
        values = [
            label,
            row.get("sample_id"),
            row.get("role"),
            row.get("country") or row.get("location") or "Unknown",
            interval,
            row.get(distance_field),
            row.get("residual"),
            row.get("reason")
            or interval_data.get("reason")
            or row.get("date_status")
            or row.get("status"),
        ]
        body.append(
            "<tr>" + "".join(f"<td>{escape(_text(value))}</td>" for value in values) + "</tr>"
        )
    headings = [
        "Sample",
        "Stable ID",
        "Role",
        "Location",
        "Collection date interval",
        units,
        "Residual (same units)",
        "Status / exclusion",
    ]
    return (
        '<div class="table-scroll" role="region" aria-label="Sample distances and dates" '
        'tabindex="0"><table><thead><tr>'
        + "".join(f'<th scope="col">{escape(label)}</th>' for label in headings)
        + "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></div>"
    )


def _panel(
    source: dict,
    diagnostic: dict,
    directory: Path,
    colours: dict[str, str],
    *,
    stem: str,
    title: str,
    anchor: str,
    description: str,
    methods: list[tuple[str, object]],
    units: str,
    slope_units: str,
) -> tuple[str, list[dict], list[tuple[str, str]]]:
    rows = _sample_rows(source, diagnostic)
    distance_field = diagnostic.get("distance_field", "distance")
    csv_name = stem + "_points.csv"
    _write_csv(directory / csv_name, rows)
    svg, png = write_date_distance_figure(
        diagnostic,
        rows,
        directory / stem,
        title=title,
        distance_label=units,
        colours=colours,
    )
    fitted = diagnostic.get("status") == "fitted"
    constant = (
        fitted and diagnostic.get("pearson_r") is None and diagnostic.get("r_squared") is None
    )
    fit_status = (
        "Constant distance" if constant else "Descriptive fit" if fitted else "Fit unavailable"
    )
    reason = (
        "Distances are constant; the slope is zero and Pearson r / R² are undefined."
        if constant
        else diagnostic.get("reason")
        or source.get("reason")
        or ("Fit available" if fitted else "No fit available")
    )
    points = _records(diagnostic.get("points"))
    caption = (
        f"{len(points)} eligible dated points; {len(rows)} retained sample rows. "
        "Horizontal bars show collection-date intervals; points use their midpoints. "
        "Bars are date precision, not confidence intervals. Residuals are observed minus fitted "
        f"distance in {units}. Location colours and input/context markers are shared across panels. "
        "Sample identities and exclusions are retained in the table and CSV."
    )
    counts = _table(
        [
            ("Retained sample rows", len(rows)),
            ("Eligible dated points", len(points)),
            ("Diagnostic status", diagnostic.get("status", source.get("status"))),
            ("Status detail", reason),
        ],
        label=title + " counts",
    )
    measures = (
        '<div class="measure-strip">'
        f"<div><small>Slope</small><b>{escape(_number(diagnostic.get('slope')))}</b>"
        f"<span>{escape(slope_units)}</span></div>"
        f"<div><small>Pearson r</small><b>{escape(_number(diagnostic.get('pearson_r'), missing='Undefined' if constant else 'Not fitted'))}</b>"
        "<span>Descriptive correlation</span></div>"
        f"<div><small>R²</small><b>{escape(_number(diagnostic.get('r_squared'), missing='Undefined' if constant else 'Not fitted'))}</b>"
        "<span>Descriptive linear fit</span></div></div>"
    )
    files = [
        (csv_name, "All sample distances, date intervals, residuals and exclusion fields"),
        (svg.name, "Scalable date-distance figure and residuals"),
        (png.name, "Raster date-distance figure and residuals"),
    ]
    panel = (
        f'<section class="stage" id="{anchor}"><div class="stage-body">'
        f'<div class="stage-head"><div><h2>{escape(title)}</h2></div>'
        f'<span class="decision review">{fit_status}</span></div>'
        f"<p>{escape(description)}</p>{counts}{measures}"
        f'<figure><a href="{svg.name}"><img src="{svg.name}" '
        f'alt="{escape(title + "; dated points, fitted line and residuals", quote=True)}"></a>'
        f'<figcaption>{escape(caption)} <a href="{svg.name}">SVG</a> · '
        f'<a href="{png.name}">PNG</a></figcaption></figure>'
        '<details class="evidence-files"><summary>Panel methods and provenance</summary>'
        + _table(methods, label=title + " methods")
        + "</details>"
        '<details class="evidence-files"><summary>Sample distances and dates</summary>'
        + _samples_table(rows, distance_field, units)
        + "</details>"
        '<details class="evidence-files"><summary>Download panel evidence</summary>'
        '<ul class="downloads">'
        + "".join(_download(*file) for file in files)
        + "</ul></details></div></section>"
    )
    return panel, [{"panel": title, **row} for row in rows], files


def write_temporal_report(
    embedding_diagnostic: dict,
    cgmlst_diagnostic: dict,
    directory: Path,
) -> Path:
    """Write the established report shell and all portable supporting artifacts."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    regression = _mapping(embedding_diagnostic.get("regression"))
    cohorts = _records(cgmlst_diagnostic.get("cohorts"))
    all_rows = _sample_rows(embedding_diagnostic, regression)
    for cohort in cohorts:
        all_rows += _sample_rows(cohort, _mapping(cohort.get("diagnostic")))
    colours = country_palette(
        sorted({str(row.get("country") or row.get("location") or "Unknown") for row in all_rows})
    )
    sections = []
    combined = []
    downloads = []
    for index, cohort in enumerate(cohorts, start=1):
        diagnostic = _mapping(cohort.get("diagnostic"))
        title = "cgMLST rooted NJ distance versus collection date"
        if len(cohorts) > 1:
            title += " · " + str(cohort.get("cohort_id", index))
        content, rows, files = _panel(
            cohort,
            diagnostic,
            directory,
            colours,
            stem=f"cgmlst_date_distance_{index}",
            title=title,
            anchor=f"cgmlst-{index}",
            description="Root-to-tip distance on the conventional cgMLST neighbour-joining tree. "
            "The root is a real sample chosen independently of collection dates.",
            methods=[
                ("Cohort", cohort.get("cohort_id")),
                ("Scheme", cohort.get("scheme_id")),
                ("Scheme version", cohort.get("scheme_version")),
                ("Database version", cohort.get("database_version")),
                ("Root sample", cohort.get("root")),
                ("Root selection", cohort.get("root_selection")),
                ("Profiles in tree", cohort.get("n_profiles")),
                ("Common callable loci", cohort.get("n_callable_loci")),
                ("Callable locus IDs", cohort.get("callable_loci")),
                ("Loci in scheme", cohort.get("scheme_loci")),
                ("Negative NJ branches", cohort.get("negative_branch_count")),
                ("Tree branch units", cohort.get("distance_units")),
                ("Fit", "Unweighted ordinary least squares using collection-date midpoints"),
            ],
            units="cgMLST mismatch fraction",
            slope_units="cgMLST mismatch fraction/year",
        )
        sections.append(content)
        combined += rows
        downloads += files
    if not cohorts:
        sections.append(
            '<section class="stage" id="cgmlst"><div class="stage-body">'
            "<h2>cgMLST rooted NJ distance versus collection date</h2>"
            '<p class="caveat">No comparable cgMLST cohort was available.</p>'
            + _table(
                [
                    ("Status", cgmlst_diagnostic.get("status")),
                    ("Exclusions", cgmlst_diagnostic.get("exclusions")),
                ],
                label="cgMLST diagnostic status",
            )
            + "</div></section>"
        )
    reference = _mapping(embedding_diagnostic.get("reference"))
    panel = _mapping(embedding_diagnostic.get("panel"))
    provenance = _mapping(embedding_diagnostic.get("provenance"))
    model = provenance.get("model", {})
    content, rows, files = _panel(
        embedding_diagnostic,
        regression,
        directory,
        colours,
        stem="protein_date_distance",
        title="Protein distance to reference versus collection date",
        anchor="protein",
        description="Mean per-locus cosine distance between each sample and one fixed real "
        "reference sample, using the same callable protein-locus panel throughout. "
        "The reference and panel are fixed independently of collection dates.",
        methods=[
            ("Reference sample", reference.get("sample_id")),
            ("Reference selection", reference.get("selection")),
            ("Reference is date independent", reference.get("date_independent")),
            ("Callable panel loci", panel.get("count", len(panel.get("loci", [])))),
            ("Panel selection", panel.get("selection")),
            ("Panel locus IDs", panel.get("loci")),
            ("Scheme", embedding_diagnostic.get("scheme")),
            ("Model / extraction provenance", model),
            ("Vector artifact", provenance.get("vector_artifact")),
            ("Embedding manifest SHA256", provenance.get("manifest_sha256")),
            ("Sample–locus mapping SHA256", provenance.get("mapping_sha256")),
            ("Fit", "Unweighted ordinary least squares using collection-date midpoints"),
        ],
        units="Mean locus cosine distance",
        slope_units="mean locus cosine distance/year",
    )
    sections.append(content)
    combined += rows
    downloads += files
    _write_csv(directory / "temporal_points.csv", combined)
    (directory / "temporal_diagnostics.json").write_text(
        json.dumps(
            {
                "embedding": embedding_diagnostic,
                "cgmlst": cgmlst_diagnostic,
                "location_colours": colours,
            },
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    counts = embedding_diagnostic.get("counts", {})
    totals = _table(
        [
            ("Prepared samples", counts.get("total_samples")),
            ("Samples with protein mappings", counts.get("mapped_samples")),
            ("Samples eligible for the fixed protein panel", counts.get("eligible_samples")),
            (
                "Eligible dated protein points",
                counts.get("dated_points", len(regression.get("points", []))),
            ),
            ("Excluded protein samples", counts.get("excluded_samples")),
        ],
        label="Prepared and eligible sample counts",
    )
    overview = (
        '<section class="stage" id="overview"><div class="stage-body">'
        "<h2>Analysis scope</h2><p>Two descriptive date-distance diagnostics use different units. "
        "They are exploratory and do not constitute a formal molecular-clock or temporal-signal test.</p>"
        + totals
        + '<p class="caveat">Protein embeddings cannot distinguish synonymous DNA changes '
        "that produce identical proteins. Embedding distances and their slopes are not substitution "
        "rates; no ancestry, ancestral date or TMRCA is inferred from these embeddings.</p></div></section>"
    )
    evidence = (
        '<section class="stage" id="evidence"><div class="stage-body"><h2>Methods and downloads</h2>'
        "<p>Regressions use eligible collection-date interval midpoints. Horizontal intervals retain "
        "the supplied date precision. Missing or invalid dates are excluded from fits; their finite "
        "distances remain in the undated panel and all sample rows remain downloadable. Too few "
        "distinct dates or constant distances are reported explicitly. The fitted line, correlation "
        "and R² are descriptive, without confidence intervals or p-values.</p>"
        '<details class="evidence-files"><summary>Exclusions</summary>'
        + _table(
            [
                ("Protein / mapping / date exclusions", embedding_diagnostic.get("exclusions", [])),
                ("cgMLST exclusions", cgmlst_diagnostic.get("exclusions", [])),
            ],
            label="Analysis exclusions",
        )
        + "</details>"
        '<details class="evidence-files" open><summary>Download all evidence</summary>'
        '<ul class="downloads">'
        + _download(
            "temporal_diagnostics.json", "Complete diagnostics, exclusions and model provenance"
        )
        + _download(
            "temporal_points.csv", "Sample rows from both diagnostics, with a named panel column"
        )
        + "".join(_download(*file) for file in downloads)
        + "</ul></details></div></section>"
    )
    cgmlst_anchor = "cgmlst-1" if cohorts else "cgmlst"
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        "<title>ChronoClade · ESM2 temporal exploration</title><style>"
        + report_styles()
        + '</style></head><body><main class="shell">'
        '<header class="identity"><div class="identity-mark">CHRONOCLADE · ESM2</div>'
        '<div class="identity-copy"><h1>Distance and collection date</h1>'
        "<p>Conventional cgMLST and experimental protein embedding diagnostics.</p>"
        '<div class="overall"><small>Analysis status</small><b>Exploratory diagnostics</b></div>'
        '</div></header><nav class="contents" aria-label="Report contents">'
        '<a href="#overview">Scope</a>'
        f'<a href="#{cgmlst_anchor}">cgMLST NJ</a><a href="#protein">Protein reference</a>'
        '<a href="#evidence">Evidence</a></nav>'
        + overview
        + "".join(sections)
        + evidence
        + '<footer class="report-close"><div>ChronoClade · Experimental ESM2 analysis · '
        "All figures and evidence are saved with this offline report.</div></footer>"
        "</main></body></html>"
    )
    output = directory / "index.html"
    output.write_text(html, encoding="utf-8")
    return output
