"""Result-only offline reports preserving the established full-report style."""

from collections import Counter
from html import escape
from pathlib import Path
import json

from chronoclade.report_components.styles import report_styles
from chronoclade.location_network.widget import interactive_network_html, widget_assets


def write_report(
    path, *, title, result, rows=(), network=None, figures=(), extra="", overview=None
):
    path = Path(path)
    counts = "".join(
        f"<h3>{escape(field.title())}</h3><p>"
        + "; ".join(
            f"{escape(str(key))}: {count}"
            for key, count in sorted(Counter(row.get(field) or "Unknown" for row in rows).items())
        )
        + "</p>"
        for field in ("country", "host")
    )
    table = (
        "<table><thead><tr><th>Genome</th><th>Role</th><th>Collection date</th><th>Country</th><th>Host</th></tr></thead><tbody>"
        + "".join(
            "<tr>"
            + f'<td title="{escape(row["sample_id"], quote=True)}">{escape(str(row.get("label") or row.get("accession") or row["sample_id"]))}</td>'
            + "".join(
                f"<td>{escape(str(row.get(field) or 'Unknown'))}</td>"
                for field in ("role", "collection_date", "country", "host")
            )
            + "</tr>"
            for row in rows
        )
        + "</tbody></table>"
    )
    plots = "".join(
        f'<figure><img src="{escape(name, quote=True)}" alt="{escape(label, quote=True)}"><figcaption>{escape(label)}</figcaption></figure>'
        for name, label in figures
        if (path.parent / name).is_file()
    )
    labels = {row["sample_id"]: row.get("label") or row["sample_id"] for row in rows}
    nearest = result.get("nearest_contexts_by_corrected_snps", [])
    if nearest:
        plots += "<h3>Nearest selected context by corrected SNP distance</h3><p>The search denominator is this exact selected assembly set. Distances are SNP counts and pair-specific callable sites.</p><table><thead><tr><th>Input</th><th>Nearest context</th><th>Corrected SNPs</th><th>Callable sites</th></tr></thead><tbody>"
        for pair in nearest:
            plots += (
                "<tr>"
                + "".join(
                    "<td>" + escape(str(value)) + "</td>"
                    for value in (
                        labels.get(pair["focal_sample"], pair["focal_sample"]),
                        labels.get(pair["context_sample"], pair["context_sample"]),
                        pair["clonal_snps"],
                        pair["callable_sites"],
                    )
                )
                + "</tr>"
            )
        plots += "</tbody></table>"
    downloads = "".join(
        f'<li><a href="{escape(p.relative_to(path.parent).as_posix(), quote=True)}">{escape(p.name)}</a></li>'
        for p in sorted(path.parent.rglob("*"))
        if p.is_file() and p != path and p.name != "job.json"
    )
    assessment = result["temporal_assessment"]
    network_html = ""
    assets = ""
    if network is not None:
        basis = "Dated-tree" if result.get("dating_status") == "dated" else "Corrected-tree"
        network_html = f"<h3>{basis} location history</h3>" + interactive_network_html(
            network, "selected-tree-network"
        )
        assets = widget_assets()
    sections = []
    if rows:
        sections.append(
            ("Selected samples", counts + '<div class="table-scroll">' + table + "</div>")
        )
    if plots or network_html:
        sections.append(("Biological tree and location history", plots + network_html))
    sections.extend(
        [
            (
                "Temporal assessment",
                f"<p><strong>{escape(assessment['code'])}</strong>: {escape(assessment['reason'])}</p><details><summary>Saved assessment</summary><pre>{escape(json.dumps(assessment, indent=2))}</pre></details>"
                + extra,
            ),
            ("Downloads", "<ul>" + downloads + "</ul>"),
        ]
    )
    if overview is None:
        overview = f"{len(result['selected_sample_ids'])} exact selected genomes"
        if result.get("branch_units"):
            overview += " · " + result["branch_units"]
        elif result.get("dating_status") == "unsupported":
            overview += " · dated tree unavailable"
    supported = assessment.get("supported", False)
    verdict_title = (
        "Dating unsupported"
        if result.get("dating_status") == "unsupported"
        else "Temporal signal supported"
        if supported
        else "Temporal signal not supported"
    )
    if assessment["code"] == "sensitivity":
        verdict_title = "Selection sensitivity assessment"
    verdict = (
        f'<div class="verdict {"supported" if supported else "not_supported"}">'
        f"<strong>{escape(verdict_title)}</strong><p>{escape(assessment['reason'])}</p></div>"
    )
    body = "".join(
        f'<section class="stage"><div class="stage-index"><span>{index}</span></div><div class="stage-body"><h2>{escape(label)}</h2>{content}</div></section>'
        for index, (label, content) in enumerate(sections, 1)
    )
    path.write_text(
        f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(title)}</title><style>{report_styles()}pre{{white-space:pre-wrap;overflow-wrap:anywhere;max-height:500px;overflow:auto}}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">CHRONOCLADE</div><div class="identity-copy"><h1>{escape(title)}</h1><p>{escape(overview)}</p>{verdict}</div></header>{body}</main>{assets}</body></html>',
        encoding="utf-8",
    )
    return path
