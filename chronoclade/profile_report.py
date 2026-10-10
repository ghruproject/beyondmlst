"""Reader-facing report for profile-first ChronoClade analyses."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any

from chronoclade.report_components.styles import report_styles
from chronoclade.sample_labels import label_analysis, read_sample_labels
from chronoclade.location_network.widget import interactive_network_html, widget_assets

from chronoclade.report import (
    _context_geography_visual,
    _context_selection_summary,
    _country_network_visual,
    _neighbourhood_visual,
    _public_health_visual,
    _recombination_visual,
)


def _profile_styles() -> str:
    """Use the established lineage-report design for profile and corrected reports."""
    return report_styles() + """
.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));border:1px solid var(--ink);margin:24px 0}.metrics>div{padding:16px;border-right:1px solid var(--line);border-bottom:1px solid var(--line);min-width:0}dt{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.06em}dd{margin:6px 0 0;font-size:21px;font-weight:700;overflow-wrap:anywhere;font-variant-numeric:tabular-nums}.muted,.empty{color:var(--muted)}.notice,.caution{background:#fff8dc;padding:15px;margin:18px 0}.table-wrap{max-width:100%;overflow:auto;scrollbar-color:var(--blue) #e5e7ec}.table-wrap table{min-width:580px}.table-wrap:focus-visible{outline:3px solid #ffb800;outline-offset:3px}.figures{display:block}.figures figure{margin:28px 0}.figures img{max-height:none;object-fit:contain}.groups{display:block}.group{padding:22px 0;border-top:1px solid var(--line)}.group h3{margin:0 0 12px}.group .metrics{margin:16px 0}.compact-group,.coverage-details{border-top:1px solid var(--line);padding:0;background:#fff}.compact-group summary,.coverage-details summary{cursor:pointer;padding:14px 0;font-weight:700}.compact-group[open]{padding-bottom:18px}.warnings{margin:8px 0}.profile-subsection{margin-top:32px}.profile-subsection h2{font-size:27px}.profile-subsection>p{margin-top:12px}@media(max-width:800px){.table-wrap table{min-width:580px}.metrics{grid-template-columns:1fr}.metrics>div{border-right:0}.table-wrap th:first-child,.table-wrap td:first-child{position:sticky;left:0;background:#fff}.contents{overflow:auto}}
    """


def _mapping(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _items(value: object) -> list[Any]:
    return value if isinstance(value, list) else []


def _text(value: object, fallback: str = "Not reported") -> str:
    if value is None or value == "":
        return fallback
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _number(value: object) -> str:
    if value is None or value == "":
        return "Not reported"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return escape(str(value))
    return f"{number:,.3g}"


def _has_reportable_value(value: object) -> bool:
    """Test emptiness without comparing NumPy scalars to containers."""
    if value is None:
        return False
    if isinstance(value, str):
        return value != ""
    if isinstance(value, (list, dict, tuple, set)):
        return bool(value)
    return True


def _records(value: object) -> list[dict[str, Any]]:
    if isinstance(value, dict):
        raw = value.get("records", value.get("items", value.get("groups", [])))
        if isinstance(raw, list):
            return [item for item in raw if isinstance(item, dict)]
        return [item for item in value.values() if isinstance(item, dict)]
    return [item for item in _items(value) if isinstance(item, dict)]


def _table(rows: list[dict[str, Any]], columns: list[tuple[str, str]], *, empty: str) -> str:
    if not rows:
        return f'<p class="empty">{escape(empty)}</p>'
    head = "".join(f"<th scope=\"col\">{escape(label)}</th>" for _, label in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{escape(_text(row.get(key)))}</td>" for key, _ in columns)
        body.append(f"<tr>{cells}</tr>")
    return (
        '<div class="table-wrap" tabindex="0" role="region" aria-label="Analysis results table">'
        f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"
    )


def _path_file(directory: Path, paths: dict[str, Any], *keys: str) -> Path | None:
    for key in keys:
        raw = paths.get(key)
        if not isinstance(raw, str) or not raw:
            continue
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = directory / candidate
        try:
            relative = candidate.resolve().relative_to(directory.resolve())
        except (OSError, ValueError):
            continue
        if (directory / relative).is_file():
            return relative
    return None


def _asset(directory: Path, paths: dict[str, Any], keys: tuple[str, ...], label: str) -> str:
    path = _path_file(directory, paths, *keys)
    if path is None:
        return f'<p class="empty">{escape(label)} is not available in this output.</p>'
    rel = path.as_posix()
    suffix = path.suffix.casefold()
    if suffix in {".svg", ".png", ".jpg", ".jpeg", ".webp"}:
        return (
            f'<figure><a href="{escape(rel, quote=True)}"><img loading="lazy" '
            f'src="{escape(rel, quote=True)}" alt="{escape(label, quote=True)}"></a>'
            f'<figcaption>{escape(label)} · <a href="{escape(rel, quote=True)}">Open full-size figure</a>'
            "</figcaption></figure>"
        )
    return f'<p><a href="{escape(rel, quote=True)}">{escape(label)} file</a></p>'


def _cohort_figures(value: object, directory: Path, field: str, label: str) -> str:
    figures = []
    for row in _records(value):
        sample_ids = row.get("sample_ids")
        if field == "pcoa_figure" and isinstance(sample_ids, list) and len(sample_ids) < 2:
            continue
        raw = row.get(field)
        if not isinstance(raw, str) or not raw:
            continue
        path = _path_file(directory, {field: raw}, field)
        if path is None:
            continue
        cohort = _text(row.get("cohort_id", "Cohort"))
        if path.suffix.casefold() not in {".svg", ".png", ".jpg", ".jpeg", ".webp"}:
            figures.append(
                f'<p><a href="{escape(path.as_posix(), quote=True)}">{escape(cohort)} · {escape(label)} file</a></p>'
            )
            continue
        figures.append(
            f'<figure><a href="{escape(path.as_posix(), quote=True)}"><img loading="lazy" '
            f'src="{escape(path.as_posix(), quote=True)}" alt="{escape(cohort + " " + label, quote=True)}"></a>'
            f"<figcaption>{escape(cohort)} · {escape(label)} · "
            f'<a href="{escape(path.as_posix(), quote=True)}">Open full size</a></figcaption></figure>'
        )
    return "".join(figures)


def _root_to_tip(value: object) -> str:
    rows = _records(value)
    if not rows:
        return '<p class="empty">An allele-unit root-to-tip screen was not available.</p>'
    content = []
    for row in rows:
        content.append(
            "<tr>"
            f"<th scope=\"row\">{escape(_text(row.get('cohort_id')))}</th>"
            f"<td>{escape(_text(row.get('root')))}</td>"
            f"<td>{escape(_number(row.get('midpoint_date_slope')))}</td>"
            f"<td>{escape(_number(row.get('midpoint_date_pearson_r')))}</td>"
            f"<td>{escape(_text(row.get('interpretation')))}</td></tr>"
        )
    return (
        '<div class="table-wrap" tabindex="0" role="region" aria-label="Exploratory root-to-tip results">'
        '<table><thead><tr><th>Cohort</th><th>Root</th><th>Slope</th><th>Pearson r</th><th>Interpretation</th></tr></thead>'
        f"<tbody>{''.join(content)}</tbody></table></div>"
    )


def _range_text(low: object, high: object) -> str:
    if low is None:
        return "Not reported"
    return str(low) if high is None or high == low else f"{low}–{high}"


def _network_tables(row: dict[str, Any], *, focus_inputs: bool) -> str:
    edges = _records(row.get("directed_edges", row.get("edges")))
    countries = set(_items(row.get("input_countries")))
    if focus_inputs:
        edges = [edge for edge in edges if _edge_countries(edge)[0] in countries or _edge_countries(edge)[1] in countries]
    edges = [{**edge,
        "country_a": _edge_countries(edge)[0], "country_b": _edge_countries(edge)[1],
        "displayed_changes": edge.get("representative_count", edge.get("displayed_changes")),
        "minimum_changes": edge.get("min_changes"), "maximum_changes": edge.get("max_changes"),
    } for edge in edges]
    directed = "directed_edges" in row
    table = _table(edges, [("country_a", "From ancestral country" if directed else "Country A"),
        ("country_b", "To descendant country" if directed else "Country B"),
        ("displayed_changes", "Changes in displayed reconstruction"),
        ("minimum_changes", "Minimum changes"), ("maximum_changes", "Maximum changes")],
        empty="No inferred country connections are available for this view.")
    nearest = []
    for edge in _records(row.get("nearest_edges", row.get("nearest_country_connections"))):
        item = dict(edge)
        item["mismatches"] = _range_text(edge.get("allele_mismatches_min"), edge.get("allele_mismatches_max"))
        item["shared_loci"] = _range_text(edge.get("shared_loci_min"), edge.get("shared_loci_max"))
        nearest.append(item)
    return (
        '<h4>Closest-relative country connections</h4>'
        '<p>These connections come from actual nearest-neighbour cgMLST comparisons, including all distance ties. Subgroup views select input genomes while retaining their closest matches across the complete selected CG pool. '
        'They do not infer movement between countries. Missing country metadata is retained as Unknown.</p>'
        + _table(nearest, [("source", "Input country"), ("target", "Nearest-relative country"),
            ("query_count", "Input genomes"), ("context_count", "Public relatives"),
            ("comparison_count", "Comparisons, including ties"),
            ("mismatches", "Allele differences"), ("shared_loci", "Jointly compared loci")],
            empty="No closest-relative country connections are available for this view.")
        + f'<details class="coverage-details"><summary>Inferred country connections ({len(edges)} pairs)</summary>'
        '<p class="muted">Ranges count changes across all equally optimal parsimony assignments on the fixed rooted NJ tree. '
        'They are exact conditional counts, not confidence or probabilities. A pair can have zero changes in the '
        'displayed representative reconstruction and still have a non-zero maximum. Ranges are calculated '
        'for each ordered pair separately; their maxima need not occur together in one reconstruction.</p>' + table + '</details>'
        + _network_metrics(row)
    )


def _network_metrics(row: dict[str, Any]) -> str:
    metrics = _records(row.get("network_metrics"))
    if not metrics:
        return ""
    metrics = [{**item, **{key: f"{item[key]:.3f}" if item.get(key) is not None else "Undefined"
                          for key in ("source_hub_ratio", "betweenness", "closeness")}}
               for item in metrics]
    return ('<details class="coverage-details"><summary>Country network metrics</summary>'
            '<p>Incoming and outgoing links count distinct country connections in the displayed history. '
            'Incoming and outgoing changes count tree branches. As in StrainHub, the source/hub ratio is outgoing links '
            'divided by incoming plus outgoing links. Centralities use unweighted links in the displayed rooted reconstruction.</p>'
            + _table(metrics, [("country", "Country"), ("in_degree", "Incoming links"),
                ("out_degree", "Outgoing links"), ("in_changes", "Incoming changes"),
                ("out_changes", "Outgoing changes"), ("degree", "Degree"),
                ("betweenness", "Betweenness"), ("closeness", "Closeness"),
                ("source_hub_ratio", "Source/hub ratio")],
                     empty="No country metrics are available.")
            + '</details>')


def _edge_countries(edge: dict[str, Any]) -> tuple[object, object]:
    return (edge.get("country_a", edge.get("source")), edge.get("country_b", edge.get("target")))


def _network_view(row: dict[str, Any], directory: Path, *, focus_inputs: bool) -> str:
    tree = ('<details class="coverage-details"><summary>Country-coloured NJ tree</summary>'
            + _figure_asset(directory, row, "country_tree_figure", "Representative country states on NJ tree")
            + '</details>')
    representative = "network_figure" if focus_inputs else "full_network_figure"
    possible = "possible_input_network_figure" if focus_inputs else "possible_network_figure"
    reconstruction = _mapping(row.get("reconstruction"))
    nodes = _records(reconstruction.get("nodes"))
    ambiguous = row.get("ambiguous_internal_node_count")
    if ambiguous is None:
        ambiguous = sum(1 for node in nodes if not node.get("is_tip")
                        and isinstance(node.get("allowed_states"), list)
                        and len(node["allowed_states"]) > 1)
    metric_values = [("Minimum total changes", reconstruction.get("optimum_changes", row.get("minimum_total_changes"))),
                     ("Ambiguous internal nodes", ambiguous),
                     ("Unknown-country samples", row.get("unknown_country_count"))]
    metrics = _metric_cards([(label, value) for label, value in metric_values if value is not None])
    roots = row.get("root_count", row.get("tested_roots"))
    roots = len(roots) if isinstance(roots, list) else roots
    root_note = ""
    if row.get("method") is not None or roots is not None:
        root_note = '<p class="muted">' + escape(_text(row.get("method")))
        if roots is not None:
            root_note += '. Separate possible-link diagnostic: ' + escape(_text(roots)) + ' roots checked'
        root_note += '.</p>'
    default_asset = _figure_asset(directory, row, representative, "Country-state transition network")
    possible_asset = _figure_asset(directory, row, possible, "Representative network with alternative possible links")
    identity = f"{row.get('cohort_id', 'cohort')}-{row.get('id', 'all')}-{'input' if focus_inputs else 'all'}"
    static = '<details class="coverage-details"><summary>Static network figure</summary>' + default_asset + '</details>'
    alternative_static = '<details class="coverage-details"><summary>Static network figure</summary>' + possible_asset + '</details>'
    return metrics + root_note + '<div data-network-representative>' + \
        interactive_network_html(row, identity + '-representative', focus_inputs=focus_inputs) + tree + static + '</div>' + \
        '<div data-network-possible hidden>' + \
        interactive_network_html(row, identity + '-possible', focus_inputs=focus_inputs, include_possible=True) + \
        tree + alternative_static + '</div>' + _network_tables(row, focus_inputs=focus_inputs)


def _figure_asset(directory: Path, row: dict[str, Any], field: str, alt: str) -> str:
    safe_row = dict(row)
    path = _path_file(directory, row, field)
    if path is None or path.suffix.casefold() not in {".svg", ".png", ".jpg", ".jpeg", ".webp"}:
        safe_row.pop(field, None)
    return _asset(directory, safe_row, (field,), alt)


def _location_network(value: object, directory: Path, cohorts: object = None) -> str:
    rows = _records(value)
    if not rows:
        return '<p class="empty">No location network reconstruction was reported.</p>'
    legacy_figures = {row.get("cohort_id"): row.get("network_figure") for row in _records(cohorts)}
    panels = []
    for index, raw in enumerate(rows):
        row = dict(raw)
        row.setdefault("network_figure", legacy_figures.get(row.get("cohort_id")))
        options = [dict(row, label="All inputs")] + _records(row.get("views"))
        rendered = [{"label": _text(view.get("label")),
            "input": _network_view(view, directory, focus_inputs=True),
            "all": _network_view(view, directory, focus_inputs=False)} for view in options]
        payload = json.dumps(rendered, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
        identifier = f"network-viewer-{index}"
        subgroup_options = "".join(f'<option value="{number}">{escape(view["label"])}</option>' for number, view in enumerate(rendered))
        panels.append(
            f'<div id="{identifier}" class="network-viewer"><h3>{escape(_text(row.get("cohort_id", "Cohort")))}</h3>'
            + '<div class="network-controls" hidden>'
            f'<label for="{identifier}-group">Input subgroup </label>'
            f'<select id="{identifier}-group" data-network-group>{subgroup_options}</select> '
            f'<label for="{identifier}-scope">Country links </label>'
            f'<select id="{identifier}-scope" data-network-scope><option value="all">All country links</option>'
            '<option value="input">Input-country links</option></select> '
            f'<label for="{identifier}-reconstruction">Reconstruction </label>'
            f'<select id="{identifier}-reconstruction" data-network-reconstruction><option value="representative">Representative reconstruction</option>'
            '<option value="possible">Include alternative possible links</option></select></div>'
            + '<p>The complete country network is shown by default. Select a country to highlight its neighbours, '
            'drag countries to separate them, or zoom and pan. Node colours show country states; dark outlines mark countries with input genomes. '
            'Arrows follow ancestral-to-descendant state changes on the displayed rooted tree. '
            'Link thickness counts these changes; dashed links mark variation across optimal assignments.</p>'
            + '<div data-network-output>' + rendered[0]["all"] + '</div>'
            + '<noscript><style>.cc-network-widget{display:none}</style><p>Open the static network figure above. Enable JavaScript for the interactive network.</p></noscript>'
            + f'<script type="application/json" data-network-views>{payload}</script></div>'
        )
    script = """<script>
(function(){document.querySelectorAll('.network-viewer').forEach(function(viewer){
var data=JSON.parse(viewer.querySelector('[data-network-views]').textContent);
var group=viewer.querySelector('[data-network-group]');var scope=viewer.querySelector('[data-network-scope]');
var reconstruction=viewer.querySelector('[data-network-reconstruction]');var output=viewer.querySelector('[data-network-output]');
function update(){var chosen=data[Number(group.value)];if(chosen&&(scope.value==='input'||scope.value==='all')){if(window.ChronoCladeNetworks)window.ChronoCladeNetworks.destroy(output);output.innerHTML=chosen[scope.value];var possible=reconstruction.value==='possible';output.querySelectorAll('div[data-network-representative]').forEach(function(el){el.hidden=possible;});output.querySelectorAll('div[data-network-possible]').forEach(function(el){el.hidden=!possible;});if(window.ChronoCladeNetworks)window.ChronoCladeNetworks.mount(output);}}
group.addEventListener('change',update);scope.addEventListener('change',update);reconstruction.addEventListener('change',update);viewer.querySelector('.network-controls').hidden=false;
});})();</script>"""
    return "".join(panels) + widget_assets() + script


def _download_links(directory: Path, paths: dict[str, Any]) -> str:
    links = []
    for label, raw in paths.items():
        if not isinstance(raw, str) or not raw:
            continue
        path = _path_file(directory, {"path": raw}, "path")
        if path is None:
            continue
        rel = path.as_posix()
        links.append(
            f'<li><a href="{escape(rel, quote=True)}">{escape(label.replace("_", " "))}</a>'
            f' <small>{escape(path.name)}</small></li>'
        )
    if not links:
        return '<p class="empty">No additional evidence files are available in this output.</p>'
    return '<ul class="downloads">' + "".join(links) + "</ul>"


def _metric_cards(values: list[tuple[str, object]]) -> str:
    return '<dl class="metrics">' + "".join(
        f"<div><dt>{escape(label)}</dt><dd>{escape(_text(value))}</dd></div>" for label, value in values
    ) + "</dl>"


def _group_cards(value: object, *, kind: str, bootstrap_requested: object = None) -> str:
    rows = _records(value)
    if not rows:
        return '<p class="empty">No groups were reported for this analysis.</p>'
    if kind in {"persistence", "concentration"}:
        return _group_details(rows, kind=kind)
    cards = []
    for i, row in enumerate(rows, 1):
        title = _text(row.get("group_id", row.get("id", row.get("label", f"Group {i}"))))
        description = row.get("interpretation", row.get("summary", row.get("description", "")))
        sample_ids = row.get("sample_ids")
        if isinstance(sample_ids, list) and "sample_count" not in row and "count" not in row:
            row = {**row, "sample_count": len(sample_ids)}
        details = []
        for key, label in (
            ("sample_count", "Samples"),
            ("count", "Samples"),
            ("first_year", "First year"),
            ("last_year", "Last year"),
            ("years", "Years represented"),
            ("observed_years", "Observed years"),
            ("observed_year_count", "Observed years"),
            ("year_span", "Year span"),
            ("records_with_time_and_place", "Records with time and place"),
            ("total_records", "Total records"),
            ("largest_cell_fraction_of_annotated", "Largest annotated cell fraction"),
            ("method", "Method"),
            ("minimum_pair_bootstrap_coassignment", "Minimum pair bootstrap coassignment"),
            ("valid_bootstrap_replicates", "Valid bootstrap replicates"),
            ("countries", "Countries"),
            ("regions", "Regions"),
            ("time_window", "Time window"),
            ("location", "Location"),
            ("concentration", "Concentration"),
            ("score", "Score"),
        ):
            value = row.get(key)
            if _has_reportable_value(value):
                details.append((label, value))
        card = f'<article class="group"><h3>{escape(title)}</h3>'
        if description:
            card += f"<p>{escape(_text(description))}</p>"
        if details:
            card += _metric_cards(details)
        cells = row.get("cells")
        if kind == "concentration" and isinstance(cells, list) and cells:
            card += _table(
                [item for item in cells if isinstance(item, dict)],
                [("location", "Location"), ("year", "Year"), ("count", "Samples")],
                empty="No time/place cells were reported.",
            )
        warning = row.get("warning")
        if warning:
            card += f'<p class="caution">{escape(_text(warning))}</p>'
        if kind == "general":
            support = row.get("minimum_pair_bootstrap_coassignment")
            valid = row.get("valid_bootstrap_replicates")
            if support is not None:
                card += f'<p class="muted">Minimum pairwise bootstrap co-assignment: {escape(_number(support))}.</p>'
            elif row.get("singleton"):
                card += '<p class="muted">Bootstrap co-assignment is not defined for a single-member group.</p>'
            elif valid == 0 and bootstrap_requested == 0:
                card += '<p class="muted">Bootstrap support was not assessed because zero replicates were requested.</p>'
            elif valid == 0:
                card += '<p class="muted">Bootstrap support is unknown because no valid replicates were available; this is not evidence of instability.</p>'
        if kind == "concentration":
            card += '<p class="muted">This describes concentration in the reported time and place; it does not show persistence across years.</p>'
        elif kind == "persistence":
            card += '<p class="muted">This describes recurrence across years; it does not imply concentration in one place.</p>'
        else:
            card += '<p class="muted">This is a descriptive grouping of the analysed profiles using the stated distance method and threshold.</p>'
        cards.append(card + "</article>")
    return '<div class="groups">' + "".join(cards) + "</div>"


def _group_details(rows: list[dict[str, Any]], *, kind: str) -> str:
    """Keep the two group summaries scannable while retaining per-group detail."""
    details = []
    for index, row in enumerate(rows, 1):
        group_id = _text(row.get("group_id", row.get("id", f"Group {index}")))
        sample_ids = row.get("sample_ids")
        sample_count = row.get("sample_count", row.get("total_records"))
        if sample_count is None and isinstance(sample_ids, list):
            sample_count = len(sample_ids)
        if kind == "persistence":
            years = row.get("observed_years", row.get("years", []))
            year_values = years if isinstance(years, list) else [years] if years else []
            year_count = row.get("observed_year_count", len(year_values))
            try:
                assessable = int(year_count) >= 2
            except (TypeError, ValueError):
                assessable = len(year_values) >= 2
            if assessable:
                summary = (
                    f"Observed across {year_count} collection years"
                    + (f" ({_text(year_values)})." if year_values else ".")
                )
                interpretation = (
                    "Observed across collection years; gaps do not establish continuous persistence. "
                    "These observations do not establish a reoccurring or persisting strain in the epidemiological sense."
                )
            elif year_values:
                summary = f"Observed in one collection year ({_text(year_values[0])}); persistence cannot be assessed."
                interpretation = "Evidence is limited to one collection year; year-to-year persistence cannot be assessed."
            else:
                summary = "No collection year was available; persistence cannot be assessed."
                interpretation = "No dated members were available to assess persistence across years."
            metrics = _metric_cards([
                ("Observed years", ", ".join(map(str, year_values)) if year_values else "None"),
                ("Year span", row.get("year_span")),
                ("Dated records", row.get("dated_records")),
                ("Total records", sample_count),
            ])
            body = f"<p>{escape(interpretation)}</p>{metrics}"
        else:
            cells = [cell for cell in _items(row.get("cells")) if isinstance(cell, dict)]
            total = row.get("records_with_time_and_place")
            largest = max(cells, key=lambda cell: cell.get("count", 0), default=None)
            if largest:
                concentration = f"Largest observed cell: {_text(largest.get('location'))}, {_text(largest.get('year'))} ({_text(largest.get('count'))} records)."
            else:
                concentration = "No time-and-place cells were available."
            summary = f"{concentration} {_text(total, 'Not reported')} records annotated."
            interpretation = _text(
                row.get("interpretation"),
                "Descriptive concentration among annotated sampled genomes; no significance or outbreak claim.",
            )
            metrics = _metric_cards([
                ("Records with time and place", total),
                ("Total records", sample_count),
                ("Largest annotated cell fraction", row.get("largest_cell_fraction_of_annotated")),
            ])
            cell_table = _table(
                cells,
                [("location", "Location"), ("year", "Year"), ("count", "Records")],
                empty="No time/place cells were reported.",
            )
            body = f"<p>{escape(interpretation)}</p>{metrics}<h4>Observed time/place cells</h4>{cell_table}"
        details.append(
            f'<details class="group compact-group"><summary><strong>{escape(group_id)}</strong> · {escape(summary)}</summary>'
            f"{body}</details>"
        )
    return '<div class="compact-groups">' + "".join(details) + "</div>"


def _cohorts(value: object, directory: Path) -> str:
    rows = _records(value)
    if not rows:
        return '<p class="empty">No cohort-level tree outputs were reported.</p>'
    rendered = []
    for row in rows:
        label = _text(row.get("cohort_id", row.get("id", row.get("name", "Cohort"))))
        units = _text(row.get("tree_units", row.get("units", "Not reported")))
        denominator = _text(row.get("coverage_denominator"))
        tree = _text(row.get("tree_path", row.get("tree", "Not reported")))
        raw_tree = row.get("tree_path", row.get("tree"))
        if isinstance(raw_tree, str):
            tree_path = _path_file(directory, {"tree": raw_tree}, "tree")
            if tree_path is not None:
                rel = tree_path.as_posix()
                tree = f'<a href="{escape(rel, quote=True)}">{escape(rel)}</a>'
        warnings = row.get("warnings", [])
        warning_html = "" if not warnings else '<ul class="warnings">' + "".join(
            f"<li>{escape(_text(item))}</li>" for item in (warnings if isinstance(warnings, list) else [warnings])
        ) + "</ul>"
        rendered.append(
            f"<tr><th scope=\"row\">{escape(label)}</th><td>{escape(units)}</td>"
            f"<td>{escape(denominator)}</td><td>{tree}</td><td>{warning_html or '—'}</td></tr>"
        )
        root = row.get("exploratory_root")
        if root:
            rendered.append(
                f'<tr><td colspan="5"><strong>Exploratory root:</strong> {escape(_text(root))}</td></tr>'
            )
    return (
        '<div class="table-wrap" tabindex="0" role="region" aria-label="Cohort tree summary">'
        '<table><thead><tr><th>Cohort</th><th>Tree units</th><th>Comparison denominator</th><th>Tree output</th><th>Warnings</th></tr></thead>'
        f"<tbody>{''.join(rendered)}</tbody></table></div>"
    )


def _coverage(analysis: dict[str, Any]) -> str:
    coverage = _mapping(analysis.get("coverage"))
    exclusions = analysis.get("exclusions", {})
    if not coverage and not exclusions:
        return '<p class="empty">Coverage and exclusions were not reported.</p>'
    rows = []
    for key, value in coverage.items():
        rows.append({"measure": key.replace("_", " "), "value": value})
    if isinstance(exclusions, list):
        excluded_counts: dict[str, int] = {}
        for item in exclusions:
            if isinstance(item, dict):
                reason = str(item.get("reason", "unspecified"))
                excluded_counts[reason] = excluded_counts.get(reason, 0) + 1
        for reason, count in sorted(excluded_counts.items()):
            rows.append({"measure": f"Excluded: {reason.replace('_', ' ')}", "value": count})
        if exclusions:
            rows.append({"measure": "Excluded sample IDs", "value": ", ".join(
                str(item.get("sample_id")) for item in exclusions if isinstance(item, dict)
            )})
    else:
        for key, value in _mapping(exclusions).items():
            rows.append({"measure": f"Excluded: {key.replace('_', ' ')}", "value": value})
    return _table(rows, [("measure", "Coverage measure"), ("value", "Count / value")], empty="Coverage unavailable")


def _input_coverage(provenance: dict[str, Any]) -> str:
    coverage = _mapping(provenance.get("coverage"))
    rows = []
    for key, label in (("queries", "Query dataset"), ("context", "Public context")):
        dataset = _mapping(coverage.get(key))
        if dataset:
            rows.append(
                {
                    "dataset": label,
                    "total": dataset.get("total"),
                    "profiles_available": dataset.get("profiles_available"),
                    "profiles_missing": dataset.get("profiles_missing"),
                }
            )
    return _table(
        rows,
        [("dataset", "Input source"), ("total", "Records"),
         ("profiles_available", "Profiles available"), ("profiles_missing", "Records without profiles")],
        empty="Separate query and public-context coverage was not supplied.",
    )


def _geography_table(value: object, *, label: str, include_origin: bool) -> str:
    rows = _records(value)
    if not rows:
        return f'<p class="empty">{escape(label)} breakdown was not available.</p>'
    totals: dict[str, int] = {}
    parsed = []
    for row in rows:
        origin = str(row.get("origin", "unspecified"))
        try:
            count = int(row.get("count", 0))
        except (TypeError, ValueError):
            count = 0
        totals[origin] = totals.get(origin, 0) + count
        parsed.append((row, origin, count))
    normalized = []
    for row, origin, count in parsed:
        normalized.append({
            "origin": {"local": "Your input genomes", "query": "Your input genomes", "focal": "Your input genomes"}.get(origin, "Public comparison genomes" if origin == "context" else origin),
            "country": row.get("country", "Unknown"),
            "region": row.get("region", "Unknown"),
            "nuts2": row.get("nuts2", "Unknown"),
            "count": count,
            "denominator": totals[origin],
        })
    columns = []
    if include_origin:
        columns.append(("origin", "Dataset"))
    columns.extend([
        ("country", "Country"), ("region", "Region"), ("nuts2", "NUTS2 region"),
        ("count", "Records"), ("denominator", "Dataset denominator"),
    ])
    scope_note = (
        "Rows are country/region/NUTS2 combinations. Counts include records with missing profiles because this table uses input metadata."
        if include_origin else
        "Rows are country/region/NUTS2 combinations from the frozen public catalogue, before the profile-analysis cap."
    )
    table = _table(normalized, columns, empty=f"{label} breakdown was not available.")
    country_count = len({str(row.get("country")) for row in rows if row.get("country") not in {None, "", "Unknown"}})
    region_count = len({str(row.get("region")) for row in rows if row.get("region") not in {None, "", "Unknown"}})
    record_count = sum(count for _, _, count in parsed)
    if not include_origin or len(rows) > 10:
        summary = (
            f"{label}: {record_count:,} records · {country_count} countries · "
            f"{region_count} regions · {len(rows)} country/region rows"
        )
        table = (
            f'<details class="coverage-details"><summary>{escape(summary)}</summary>{table}'
            f'<p class="muted">{escape(scope_note)}</p></details>'
        )
    else:
        table += f'<p class="muted">{escape(scope_note)}</p>'
    return table


def _nearest_neighbour_table(rows: list[dict[str, Any]]) -> str:
    normalized = []
    for row in rows:
        item = dict(row)
        item["distance"] = _number(row.get("distance"))
        item["cohort_id"] = row.get("cohort_id") or row.get("cohort")
        try:
            item["allele_comparison"] = (
                f"{int(row.get('allele_differences')):,} / {int(row.get('shared_called_loci')):,}"
            )
        except (TypeError, ValueError):
            item["allele_comparison"] = "Not reported"
        normalized.append(item)
    return _table(
        normalized,
        [("query_id", "Query"), ("context_id", "Nearest relative"),
         ("distance", "Normalized distance"), ("allele_comparison", "Allele differences / shared loci"),
         ("call_overlap", "Call overlap"), ("cohort_id", "Comparison group"), ("country", "Country"), ("year", "Year")],
        empty="No nearest relatives were reported.",
    )


def _shared_tree_selection_summary(data: dict[str, Any]) -> str:
    selection = _mapping(data.get("shared_selection"))
    if not selection:
        return ""
    queries = len(_items(selection.get("query_ids")))
    contexts = len(_items(selection.get("selected_context_ids")))
    text = (f'<p>The displayed NJ tree, full analysis and finish analysis use one shared selection: '
            f'<strong>{queries} input genomes and {contexts} public context genomes</strong>. '
            'Nearest-neighbour comparisons, PCoA and the country network still use the complete usable profile pool.</p>')
    missing = _items(_mapping(data.get("tree_display_selection")).get("without_comparable_profiles"))
    if missing:
        text += (f'<p class="caution">{len(missing)} selected genomes have no comparable cgMLST profiles '
                 'and cannot appear in the NJ view. They remain in the shared selection for assembly analysis.</p>')
    return text


def _adaptive_context_summary(provenance: dict[str, Any]) -> str:
    audit = _mapping(provenance.get("adaptive_context_selection"))
    datasets = _records(audit.get("datasets"))
    if not datasets:
        return ""
    overview = _table(datasets, [
        ("mlst_st", "MLST ST"), ("clonal_group", "Clonal group"),
        ("input_count", "Your input genomes"),
        ("public_cg_pool_count", "Available public CG comparisons"),
        ("selected_context_count", "Chosen public context pool"),
    ], empty="No typed input datasets were available.")
    rows = []
    for dataset in datasets:
        for subgroup in _records(dataset.get("subgroups")):
            counts = _mapping(subgroup.get("context_counts"))
            prefix = subgroup.get("input_prefix", [])
            rows.append({
                "cg": dataset.get("clonal_group"),
                "prefix": ",".join(str(part) for part in prefix) if isinstance(prefix, list) else prefix,
                "inputs": subgroup.get("input_count"),
                "level5": counts.get("5"), "level6": counts.get("6"), "level7": counts.get("7"),
                "selected_level": subgroup.get("selected_level"),
                "selected_count": subgroup.get("selected_context_count"),
                "coverage": "Limited context" if subgroup.get("limited_context") else "Minimum met",
            })
    selection = _table(rows, [
        ("cg", "Clonal group"), ("prefix", "Input level-5 prefix"), ("inputs", "Inputs"),
        ("level5", "Level 5 context"), ("level6", "Level 6 context"),
        ("level7", "Level 7 context"), ("selected_level", "Chosen level"),
        ("selected_count", "Chosen context"), ("coverage", "Context coverage"),
    ], empty="No subgroup selection counts were supplied.")
    minimum = escape(_text(audit.get("min_context")))
    return (
        '<h3>Input datasets and their public context</h3>'
        '<p>Inputs are separated by MLST sequence type, then by clonal group. '
        'Each dataset is analysed separately.</p>' + overview
        + '<h3>Which closer LIN groups were chosen?</h3>'
        + '<p>For each input level-5 subgroup, start at level 7 and widen to level 6 '
        'or 5 until every input has at least ' + minimum + ' public comparisons. '
        'This minimum is a sampling setting, not a biological cutoff. Sparse groups '
        'remain visible as limited context. Counts at each level combine the matching '
        'prefixes represented by the inputs; overlapping comparisons are counted once '
        'in the chosen dataset pool.</p>' + selection
        + '<p>LIN prefixes define where to search. Closest relatives are ranked by '
        'actual cgMLST allele differences, with jointly compared loci and distance ties '
        'reported below. One set of figures describes each dataset; separate figures '
        'are not generated for every LIN level.</p>'
    )


def _adaptive_context_metadata(provenance: dict[str, Any]) -> str:
    datasets = _records(_mapping(provenance.get("adaptive_context_selection")).get("datasets"))
    blocks = []
    for dataset in datasets:
        countries = _records(dataset.get("selected_context_geography"))
        years = _records(dataset.get("selected_context_years"))
        if not countries and not years:
            continue
        blocks.append('<h3>Complete chosen context pool: '
                      + escape(_text(dataset.get("analysis_dataset"))) + '</h3>')
        blocks.append('<p>These counts use all public metadata in the chosen LIN groups, '
                      'before tree subsampling or profile exclusions. Missing locations and '
                      'collection years remain in the denominator.</p>')
        blocks.append(_table(countries, [("country", "Country"), ("region", "Region"),
                      ("nuts2", "NUTS2 region"), ("count", "Public samples"),
                      ("denominator", "Chosen context pool")],
                      empty="Country metadata for the complete chosen pool was not available."))
        blocks.append(_table(years, [("year", "Collection year"), ("count", "Public samples"),
                      ("denominator", "Chosen context pool")],
                      empty="Collection-year metadata for the complete chosen pool was not available."))
    return "".join(blocks)


def _context_funnel(provenance: dict[str, Any]) -> str:
    rows = []
    catalogue = _mapping(provenance.get("catalogue_deduplication"))
    context = _mapping(provenance.get("context_deduplication"))
    if catalogue:
        rows.extend([
            {"stage": "Public catalogue records before deduplication", "count": catalogue.get("raw_records")},
            {"stage": "Public catalogue sample units after deduplication", "count": catalogue.get("retained_units")},
        ])
    if provenance.get("catalogue_metadata_count") is not None:
        rows.append({"stage": "Frozen catalogue metadata units", "count": provenance["catalogue_metadata_count"]})
    if context:
        rows.extend([
            {"stage": "Eligible public context after query/focal exclusions and deduplication", "count": context.get("retained_units")},
        ])
    for key, label in (
        ("eligible_public_context_count", "Eligible public context records"),
        ("bounded_public_context_count", "Public context records in profile analysis"),
        ("profile_limit", "Configured public profile limit"),
    ):
        if provenance.get(key) is not None:
            rows.append({"stage": label, "count": provenance[key]})
    selection = _mapping(provenance.get("context_pool_selection"))
    audits = [_mapping(_mapping(item).get("audit")) for item in _items(selection.get("lineages"))]
    for audit in audits:
        counts = _mapping(audit.get("counts"))
        for key, label in (("priority_selected", "Selected using matching lineage-group evidence"),
                           ("background_selected", "Selected as same-ST background comparisons")):
            if key in counts:
                rows.append({"stage": label, "count": counts[key]})
    result = _table(rows, [("stage", "Context funnel step"), ("count", "Records / units")],
                    empty="Public catalogue selection counts were not supplied.")
    if any("no_comparable_evidence_same_st_fallback" in _mapping(a.get("query_status")).values()
           for a in audits):
        result += (
            '<p>Matching lineage-group evidence was unavailable for some inputs when the public pool '
            'was selected. Those comparisons were sampled from the same MLST sequence type across '
            'countries and years. They were not selected by matching cgLIN or HierCC groups.</p>'
        )
    return result


def _provenance_notes(provenance: dict[str, Any]) -> str:
    notes = []
    if provenance.get("analysis_export_status") == "unavailable_no_credentials":
        notes.append(
            "Live Pathogenwatch typing exports were unavailable because no API credentials were configured."
        )
    exports = _items(provenance.get("analysis_exports"))
    for item in exports:
        if isinstance(item, dict) and item.get("status") == "unavailable":
            notes.append(
                f"{_text(item.get('download'), 'A public typing export')} was unavailable: "
                f"{_text(item.get('reason'), 'reason not reported')}"
            )
    skipped = _items(provenance.get("context_discovery_skipped"))
    for item in skipped:
        if isinstance(item, dict):
            notes.append(
                f"Context discovery for {_text(item.get('species'))} {_text(item.get('lineage'))} "
                f"was skipped: {_text(item.get('reason'))}."
            )
    if not notes:
        return ""
    return '<aside class="notice"><strong>Typing and context availability</strong><ul>' + "".join(
        f"<li>{escape(note)}</li>" for note in notes
    ) + "</ul></aside>"


def _available_neighbours(
    value: object,
    directory: Path,
    paths: dict[str, Any],
    *,
    public_context_count: object = None,
) -> str:
    if isinstance(value, (str, Path)):
        path = _path_file(directory, {"neighbours": str(value)}, "neighbours")
        if path is None:
            return '<p class="empty">Nearest-relative results are not available in this output.</p>'
        if path.suffix.casefold() == ".json":
            try:
                raw = json.loads((directory / path).read_text(encoding="utf-8"))
                value = raw.get("records", raw.get("nearest_neighbours", raw)) if isinstance(raw, dict) else raw
            except (OSError, json.JSONDecodeError):
                return '<p class="empty">Nearest-relative file could not be read.</p>'
        else:
            return f'<p><a href="{escape(path.as_posix(), quote=True)}">Download nearest-relative results</a></p>'
    rows = _records(value)
    if not rows:
        found = _path_file(directory, paths, "nearest_neighbours", "nearest_neighbours_tsv", "nearest_neighbours_json")
        if found:
            return f'<p><a href="{escape(found.as_posix(), quote=True)}">Download nearest-relative results</a></p>'
        return '<p class="empty">Nearest-relative results were not generated for this stage.</p>'
    matched = [row for row in rows if row.get("status") == "matched" or row.get("context_id")]
    unmatched = [row for row in rows if row not in matched]
    if not matched:
        query_count = len({row.get("query_id") for row in unmatched if row.get("query_id")})
        try:
            no_context = int(public_context_count) == 0
        except (TypeError, ValueError):
            no_context = False
        empty_message = (
            "No public context profiles were included in this analysis."
            if no_context
            else "No comparable public context profiles were available for these queries."
        )
        result = (
            f'<p class="empty">{escape(empty_message)} '
            "No nearest-relative ranking can be reported.</p>"
        )
        if unmatched:
            result += (
                f'<details class="coverage-details"><summary>Query-level status ({query_count} queries)</summary>'
                + _table(unmatched, [("query_id", "Query"), ("status", "Status")],
                         empty="No query-level statuses were reported.")
                + "</details>"
            )
        return result
    result = _nearest_neighbour_table(matched)
    if unmatched:
        query_count = len({row.get("query_id") for row in unmatched if row.get("query_id")})
        result += (
            f'<details class="coverage-details"><summary>Queries without a comparable public context profile ({query_count})</summary>'
            + _table(unmatched, [("query_id", "Query"), ("status", "Status")],
                     empty="No query-level statuses were reported.")
            + "</details>"
        )
    return result


def write_profile_report(
    analysis: dict,
    *,
    directory: Path,
    stage: str = "fast",
    provenance: dict | None = None,
    report_label: str = "Fast profile report",
) -> Path:
    """Write an offline HTML report, with the profile summary before relatives."""

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    data = label_analysis(_mapping(analysis), _mapping(analysis.get("sample_labels")))
    prov = _mapping(provenance)
    paths = _mapping(data.get("paths"))
    focal = data.get("records_count", data.get("sample_count", "Not reported"))
    heading = _text(data.get("title", "Genome comparison"))
    species = data.get("species")
    if species:
        heading = f"{_text(species).replace('_', ' ')} {_text(data.get('lineage'), '')}".strip()
    identity_heading = (
        f"<i>{escape(_text(species).replace('_', ' '))}</i><span>{escape(_text(data.get('lineage'), ''))}</span>"
        if species else escape(heading)
    )
    message = ""
    raw_profiles = _mapping(data.get("coverage")).get(
        "available_profiles", _mapping(data.get("coverage")).get("profiles_available")
    )
    try:
        input_count = int(focal or 0)
        profile_count = int(raw_profiles)
        empty_input = input_count == 0 or profile_count == 0
        partial_profiles = profile_count < input_count
    except (TypeError, ValueError):
        empty_input = False
        partial_profiles = not data.get("profiles_available", True)
    if empty_input:
        message = '<aside class="notice"><strong>No usable profiles were available.</strong> '
        message += "The report retains input and coverage information so missing typing is visible; no group or neighbour findings can be inferred.</aside>"
    elif partial_profiles:
        message = '<aside class="notice"><strong>Profile coverage is incomplete.</strong> '
        message += "Interpret group and neighbour results against the coverage and exclusion counts below.</aside>"
    warning_values = _items(data.get("warnings"))
    warnings = "" if not warning_values else '<aside class="notice"><strong>Analysis notes</strong><ul>' + "".join(
        f"<li>{escape(_text(item))}</li>" for item in warning_values
    ) + "</ul></aside>"

    query_coverage = _mapping(_mapping(prov.get("coverage")).get("queries"))
    query_count = query_coverage.get("total")
    context_coverage = _mapping(_mapping(prov.get("coverage")).get("context"))
    public_context_count = context_coverage.get("total")
    public_profiles_available = context_coverage.get("profiles_available")
    context_description = prov.get("public_typing") or data.get("public_typing")
    if not context_description and public_context_count is not None:
        context_description = (
            f"{public_profiles_available} profile{'s' if public_profiles_available != 1 else ''} available / "
            f"{public_context_count} records"
            if public_profiles_available is not None
            else f"{public_context_count} records included"
        )
    if not context_description:
        context_description = "Frozen local file supplied" if "public_typing" in _mapping(prov.get("inputs")) else "Not reported"
    summary_values = [
        ("Your input genomes", query_count if query_count is not None else "Not reported"),
        ("Public comparison genomes", public_context_count if public_context_count is not None else "Not reported"),
        ("Usable profiles across both sets", raw_profiles if raw_profiles is not None else "Not reported"),
    ]
    input_description = (
        f"This report investigates your {query_count} input genomes and compares them with "
        f"{public_context_count} additional public genomes. The selection counts and methods "
        "are recorded below."
        if query_count is not None and public_context_count is not None else
        "Your input genomes are the samples being investigated. Public comparison genomes "
        "are additional records selected from the source database to provide context."
    )
    hero_scope = (
        f"{query_count} input genomes · {public_context_count} public comparisons · {report_label}"
        if query_count is not None and public_context_count is not None else
        f"{focal} genomes · Profile comparison"
    )
    adaptive_selection = bool(_records(_mapping(prov.get("adaptive_context_selection")).get("datasets")))
    comparison_scope = (
        "All public genomes in the chosen LIN context groups, before tree subsampling."
        if adaptive_selection else "Additional genomes selected for comparison."
    )
    catalogue_heading = "Full public clonal-group pool" if adaptive_selection else "Frozen public catalogue"
    figure_scope = (
        "Country and collection-year tables describe the complete chosen context pool, including records without usable profiles. "
        "PCoA and nearest-neighbour comparisons use usable profiles; the displayed NJ tree uses the shared full/finish genome selection. "
        "The full public clonal-group pool has a separate denominator. Unknown countries and regions remain visible."
        if adaptive_selection else
        "Country figures show profile-available members of each complete-comparability cohort and combine query and public context records. "
        "The tables above separately summarize resolved input metadata by origin, including rows without profiles, and the full frozen public catalogue. "
        "These pools have different denominators. Unknown country and region values remain visible."
    )
    metadata_geography = _records(data.get("metadata_geography"))
    input_geography = [row for row in metadata_geography if row.get("origin") in {"local", "query", "focal"}]
    comparison_geography = [row for row in metadata_geography if row not in input_geography]
    geography = _mapping(data.get("geography", data.get("country_breakdown")))
    geography_rows = _records(geography) or _records(data.get("country_breakdown"))
    nearest = data.get("nearest_neighbours", data.get("nearest_relative"))
    groups = data.get("genetic_groups")
    sensitivity_rows = [dict(cohort_id=cohort.get("cohort_id"), **row)
                        for cohort in data.get("cohorts", [])
                        for row in cohort.get("group_sensitivity", {}).get("thresholds", [])]
    sensitivity_table = _table(sensitivity_rows,
        [("cohort_id", "Comparison cohort"), ("distance_fraction", "Distance cutoff"),
         ("group_count", "Groups"), ("singletons", "Single genomes"),
         ("largest_group", "Largest group"), ("selected", "Selected cutoff")],
        empty="Threshold sensitivity is unavailable.") if sensitivity_rows else ""
    concentration = data.get("time_place_concentration")
    root_tip = data.get("root_to_tip")
    root_tip_text = "Allele-distance root-to-tip is exploratory and does not pass the final clock gate."
    if isinstance(root_tip, dict):
        root_tip_text += " " + _text(root_tip.get("interpretation", root_tip.get("summary", "")), "")
    from chronoclade.report_components.ordination import ordination_viewer
    pcoa_figures = ordination_viewer(data.get("cohorts"), directory)
    if not pcoa_figures:
        pcoa_figures = _cohort_figures(data.get("cohorts"), directory, "pcoa_figure", "Profile PCoA")
    if not pcoa_figures:
        pcoa_figures = _asset(directory, paths, ("pcoa", "pcoa_svg", "pcoa_png"), "Profile PCoA")
    tree_figures = _cohort_figures(data.get("cohorts"), directory, "tree_figure", "Neighbour-joining tree")
    if not tree_figures:
        tree_figures = _asset(directory, paths, ("nj_tree", "nj_tree_svg", "nj_tree_png", "tree"), "Neighbour-joining tree")
    country_figures = _cohort_figures(data.get("cohorts"), directory, "country_figure", "Country composition")
    if not country_figures:
        country_figures = '<p class="empty">Country and region figures are not available in this output.</p>'
    legacy_geography = (
        _table(
            geography_rows,
            [("country", "Country"), ("region", "Region"), ("count", "Samples"),
             ("percentage", "Percent"), ("denominator", "Denominator")],
            empty="No country or region table was available.",
        )
        if geography_rows
        else ""
    )
    sections = [
        '<main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · ANALYSIS REPORT</div>'
        f'<div class="identity-copy"><h1>{identity_heading}</h1><p>{escape(hero_scope)}</p></div></header>'
        '<nav class="contents" aria-label="Report topics"><a href="#summary">Overview</a><a href="#geography">Countries</a><a href="#nearest">Closest relatives</a><a href="#network">Country network</a><a href="#groups">Groups</a><a href="#root-to-tip">Dates</a><a href="#downloads">Methods &amp; files</a></nav>',
        f'<section class="stage" id="summary"><div class="stage-body"><h2>Your results at a glance</h2>{message}{warnings}{_metric_cards(summary_values)}'
        f'<p>{escape(input_description)}</p><p>Public comparison typing: {escape(context_description)}.</p>'
        f'{_adaptive_context_summary(prov)}'
        f'{_shared_tree_selection_summary(data)}'
        '<p>Start with the countries and closest relatives below. The comparison uses differences in shared core genes (cgMLST). Genetic relationships and concentration in time and place are reported separately.</p></div></section>',
        f'<section class="stage" id="geography"><div class="stage-body"><h2>Collection locations</h2>'
        f'{_adaptive_context_metadata(prov)}'
        f'<h3>Your input genomes</h3><p>The samples you supplied for investigation.</p>{_geography_table(input_geography, label="Your input genomes", include_origin=True)}'
        f'<h3>Public comparison genomes</h3><p>{escape(comparison_scope)}</p>{_geography_table(comparison_geography, label="Public comparison genomes", include_origin=True)}'
        f'<h3>{escape(catalogue_heading)}</h3>{_geography_table(data.get("public_catalogue_geography"), label=catalogue_heading, include_origin=False)}'
        f'{legacy_geography}'
        f'{country_figures}'
        f'<p class="muted">{escape(figure_scope)}</p></div></section>',
        f'<section class="stage" id="nearest"><div class="stage-body"><h2>Closest relatives</h2><p>These rankings use the shared cgMLST loci in the analysed profiles.</p>{_available_neighbours(nearest, directory, paths, public_context_count=public_context_count)}</div></section>',
        '<section class="stage" id="figures"><div class="stage-body"><h2>Genetic relationships</h2><div class="figures">'
        + pcoa_figures
        + tree_figures
        + "</div><p class=\"muted\">The neighbour-joining tree is a profile-distance view. Branches do not represent time or prove transmission.</p></div></section>",
        f'<section class="stage" id="network"><div class="stage-body"><h2>Country connections</h2><p>Country reconstruction depends on the rooted tree and supplied metadata. This network describes the selected LIN context pool. It uses all usable profiles, before tree-display subsampling. Input subgroups can have different context breadth, so the network does not describe the distribution of the complete ST.</p>{_location_network(data.get("location_network"), directory, data.get("cohorts"))}</div></section>',
        f'<section class="stage" id="groups"><div class="stage-body"><h2>Genetic groups</h2><p>{escape(_text(_mapping(groups).get("interpretation", "Groups are descriptive summaries of the reported profile distances.")))}</p>'
        f'<p class="muted">Groups use the explicit distance cutoff. Threshold sensitivity and locus resampling describe grouping consistency; thresholds are not biologically validated.</p>{sensitivity_table}'
        '<details class="evidence-files"><summary>Inspect genetic groups and bootstrap results</summary>'
        f'{_group_cards(groups, kind="general", bootstrap_requested=data.get("bootstrap_replicates"))}</details></div></section>',
        f'<section class="stage" id="concentration"><div class="stage-body"><h2>Concentration in time and place</h2><p>Groups concentrated in a particular time window and location are described here. These observations do not establish an outbreak.</p>{_group_cards(concentration, kind="concentration")}</div></section>',
        f'<section class="stage" id="root-to-tip"><div class="stage-body"><h2>Exploratory root-to-tip screen</h2><p>{escape(root_tip_text.strip())}</p>{_root_to_tip(root_tip)}</div></section>',
        f'<section class="stage" id="coverage"><div class="stage-body"><h2>Coverage and exclusions</h2>{_provenance_notes(prov)}<h3>Query and context inputs</h3>{_input_coverage(prov)}'
        f'<h3>Public context selection</h3>{_context_funnel(prov)}'
        f'<h3>Profile analysis</h3>{_coverage(data)}</div></section>',
        f'<section class="stage" id="cohorts"><div class="stage-body"><h2>Cohort trees</h2>{_cohorts(data.get("cohorts"), directory)}</div></section>',
        f'<section class="stage" id="downloads"><div class="stage-body"><h2>Available evidence files</h2>{_download_links(directory, paths)}</div></section>',
        '<footer class="report-close"><div>Generated by ChronoClade. Profile distances, group summaries and source tables are preserved in the supporting files.</div></footer></main>',
    ]
    document = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><meta name="color-scheme" content="light">
<title>""" + escape(heading) + """</title><style>
""" + _profile_styles() + """
</style></head><body>""" + "\n".join(sections) + "</body></html>"
    output = directory / "profile_report.html"
    output.write_text(document, encoding="utf-8")
    return output


def write_fast_group_index(
    directory: Path,
    datasets: list[dict[str, Any]],
    *,
    provenance: dict | None = None,
) -> Path:
    """Link separate input-group analyses without duplicating their figures."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    audit = _mapping(_mapping(provenance).get("adaptive_context_selection"))
    audit_datasets = _records(audit.get("datasets"))
    rows = []
    for dataset in datasets:
        lineage = _text(dataset.get("lineage"))
        matching = next((item for item in audit_datasets
                         if item.get("analysis_dataset") == lineage), {})
        candidate = dataset.get("report")
        link = '<span class="muted">Report unavailable</span>'
        if isinstance(candidate, (str, Path)):
            path = Path(candidate)
            if not path.is_absolute():
                path = directory / path
            try:
                relative = path.resolve().relative_to(directory.resolve())
            except (OSError, ValueError):
                relative = None
            if relative is not None and path.is_file():
                link = f'<a href="{escape(relative.as_posix(), quote=True)}">Open group report</a>'
        values = [dataset.get("species"), lineage,
                  matching.get("mlst_st", dataset.get("mlst_st")),
                  matching.get("clonal_group", dataset.get("clonal_group")),
                  dataset.get("input_count", matching.get("input_count")),
                  matching.get("public_cg_pool_count"),
                  dataset.get("context_count", matching.get("selected_context_count"))]
        cells = "".join(f"<td>{escape(_text(value))}</td>" for value in values)
        rows.append(f"<tr>{cells}<td>{link}</td></tr>")
    headers = ["Species", "Dataset", "MLST ST", "Clonal group", "Your inputs",
               "Available public CG pool", "Chosen public context", "Results"]
    table = ('<div class="table-wrap" tabindex="0" role="region" aria-label="Input analysis datasets">'
             '<table><thead><tr>'
             + "".join(f'<th scope="col">{header}</th>' for header in headers)
             + '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>')
    if not rows:
        table = '<p class="empty">No group reports are available yet.</p>'
    document = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>ChronoClade fast group reports</title><style>' + _profile_styles() + '</style></head>'
        '<body><main class="shell"><header class="identity">'
        '<div class="identity-mark">ChronoClade · ANALYSIS REPORT</div>'
        '<div class="identity-copy"><h1>Fast profile results</h1>'
        '<p>Choose an input dataset to explore its context.</p></div></header>'
        '<section class="stage" id="summary"><div class="stage-body">'
        '<h2>Your input datasets</h2><p>Inputs are separated by MLST sequence type, '
        'then by clonal group. Each group has its own context selection, nearest-neighbour '
        'results, PCoA and genetic tree.</p>' + table
        + '<p>Available public CG pool counts exclude your inputs. Chosen public context '
        'counts describe the complete selected LIN groups, before tree subsampling. '
        'Open a group report to see which LIN levels were chosen and the country and '
        'collection-year breakdowns.</p></div></section></main></body></html>'
    )
    output = directory / "profile_report.html"
    output.write_text(document, encoding="utf-8")
    return output


def write_stage_index(output: Path, stages: dict) -> Path:
    """Write a compact HTML index linking only to stage reports that exist."""

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    entries = []
    for name, raw in stages.items():
        if isinstance(raw, dict):
            candidate = raw.get("report", raw.get("path"))
            status = raw.get("status", "available")
        else:
            candidate, status = raw, "available"
        if not isinstance(candidate, (str, Path)):
            continue
        path = Path(candidate)
        if not path.is_absolute():
            path = output / path
        try:
            rel = path.resolve().relative_to(output.resolve())
        except (OSError, ValueError):
            continue
        if not path.is_file():
            continue
        entries.append(
            f'<li><a href="{escape(rel.as_posix(), quote=True)}">{escape(str(name))}</a>'
            f' <span>{escape(str(status))}</span></li>'
        )
    links = "".join(entries) or '<li class="empty">No stage reports are available yet.</li>'
    page = (
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>ChronoClade analysis reports</title><style>' + _profile_styles() + '</style><body><main class="shell">'
        '<header class="identity"><div class="identity-mark">ChronoClade · RUN SUMMARY</div><div class="identity-copy"><h1>Analysis reports</h1><p>Choose a completed analysis.</p></div></header>'
        '<section class="stage"><div class="stage-body"><ul>' + links + '</ul></div></section></main></body></html>'
    )
    index = output / "index.html"
    index.write_text(page, encoding="utf-8")
    return index


def write_corrected_report(report: dict, *, directory: Path) -> Path:
    """Write the corrected-tree report before any date-based analysis is run."""

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    data = label_analysis(_mapping(report), read_sample_labels(directory))
    species = _text(data.get("species"))
    lineage = _text(data.get("lineage"))
    context = _mapping(data.get("context"))
    count = data.get("sample_count", "Not reported")
    selection = _context_selection_summary(directory)
    geography = _context_geography_visual(directory, species)
    recombination = _recombination_visual(data, directory)
    neighbours = _neighbourhood_visual(directory)
    network = _country_network_visual(directory)
    public_summary, public_evidence = _public_health_visual(data, directory)

    temporal_status = _text(data.get("temporal_status", "not_assessed"))
    temporal_reason = _text(
        data.get(
            "temporal_reason",
            data.get(
                "temporal_signal_reason",
                data.get("reason", "Dating is assessed only in finish stage"),
            ),
        )
    )
    if temporal_status == "not_assessed":
        temporal_section = (
            '<section class="stage" id="dating"><div class="stage-body"><h2>Dating</h2><p><strong>Not assessed in this stage.</strong> '
            f"{escape(temporal_reason)}</p><p>No temporal test or calendar tree is implied by this report.</p></div></section>"
        )
    else:
        temporal_section = (
            '<section class="stage" id="dating"><div class="stage-body"><h2>Dating status</h2>'
            f"<p>{escape(temporal_status)}: {escape(temporal_reason)}</p></div></section>"
        )

    context_values = [
        ("Analysed genomes", count),
        ("Focal samples", context.get("local_samples", context.get("focal_samples", "Not reported"))),
        ("Selected context", context.get("context_samples", context.get("selected_context", "Not reported"))),
        ("Context locations", context.get("context_locations", "Not reported")),
    ]
    files = _mapping(data.get("outputs"))
    output_links = []
    for label, key in (
        ("Corrected genetic tree", "genetic_tree"),
        ("Pairwise SNP distances", "clonal_pairwise_distances"),
        ("Nearest-relative results", "nearest_neighbours"),
        ("Country network data", "country_network"),
    ):
        raw = files.get(key)
        if not isinstance(raw, str):
            continue
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = directory / candidate
        try:
            relative = candidate.resolve().relative_to(directory.resolve())
        except (OSError, ValueError):
            continue
        if candidate.is_file():
            output_links.append(
                f'<li><a href="{escape(relative.as_posix(), quote=True)}">{escape(label)}</a></li>'
            )
    stage_link = ""
    fast_candidates = [directory / prefix / name for prefix in ("../../fast", "../../../../fast") for name in ("profile_report.html", "index.html")]
    for candidate in fast_candidates:
        if candidate.is_file():
            from os.path import relpath
            relative = Path(relpath(candidate.resolve(), directory.resolve()))
            stage_link = (
                f'<p><a href="{escape(relative.as_posix(), quote=True)}">Open the earlier profile-first report</a></p>'
            )
            break
    links_section = (
        "<ul>" + "".join(output_links) + "</ul>" if output_links else
        '<p class="muted">No additional output files were listed in the analysis record.</p>'
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><meta name="color-scheme" content="light">
<title>ChronoClade corrected report — {escape(species)} {escape(lineage)}</title><style>
{_profile_styles()}
</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · ANALYSIS REPORT</div><div class="identity-copy">
<h1><i>{escape(species.replace("_", " "))}</i><span>{escape(lineage)}</span></h1><p>{escape(_text(count))} genomes · Recombination-adjusted analysis</p>{stage_link}</div></header>
<nav class="contents" aria-label="Report topics"><a href="#summary">Overview</a><a href="#countries">Countries</a><a href="#relationships">Closest relatives</a><a href="#network">Country network</a><a href="#interpretation">Interpretation</a><a href="#dating">Dating</a><a href="#outputs">Methods &amp; files</a></nav>
<section class="stage" id="summary"><div class="stage-body"><h2>Your results at a glance</h2>{_metric_cards(context_values)}<p class="muted">Country composition and corrected-tree counts use the records available to this analysis. Public-catalogue totals, where reported, have a separate denominator.</p></div></section>
<section class="stage" id="countries"><div class="stage-body"><h2>Collection locations</h2>{selection}{geography}<p class="muted">Recorded locations describe submitted metadata. They do not establish where infection was acquired or population prevalence.</p></div></section>
<section class="stage" id="relationships"><div class="stage-body"><h2>Closest relatives</h2>{neighbours}</div></section>
<section class="stage" id="recombination"><div class="stage-body"><h2>Recombination evidence</h2>{recombination}</div></section>
<section class="stage" id="interpretation"><div class="stage-body"><h2>Genomic interpretation</h2>{public_summary}{public_evidence}</div></section>
<section class="stage" id="network"><div class="stage-body"><h2>What location changes does the tree suggest?</h2>{network}<p class="muted">Location changes depend on the rooted tree and supplied sample metadata. They include reconstruction uncertainty and are not proof of transmission or acquisition direction.</p></div></section>
{temporal_section}
<section class="stage" id="outputs"><div class="stage-body"><h2>Available files</h2>{links_section}</div></section><footer class="report-close"><div>Generated by ChronoClade. The supporting files preserve the corrected tree and comparison results.</div></footer></main></body></html>"""
    output = directory / "corrected_report.html"
    output.write_text(document, encoding="utf-8")
    return output
