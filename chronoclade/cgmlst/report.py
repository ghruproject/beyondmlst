"""The existing report shell presents saved profile partitions and selections."""

import html
from chronoclade.report_components.styles import report_styles


def write_partition_index(
    products, audit, *, records_count, lin_level, hiercc_level, status, output
):
    items = "".join(
        f"<tr><td>{html.escape(row['label'])}</td><td>{row['input_count']}</td>"
        f'<td>{row["context_count"]}</td><td><a href="{html.escape(row["report"])}">Report</a> · '
        f'<a href="{html.escape(row["ensemble"])}">Selections</a></td></tr>'
        for row in products
    )
    unresolved = "".join(
        f"<li>{html.escape(row['sample_id'])}: {row['reason']}</li>"
        for row in audit["unresolved_inputs"]
    )
    counts = "".join(
        f"<tr><td>{row['level']}</td><td>ST{html.escape(row['mlst_st'])}</td>"
        f"<td>{html.escape('.'.join(map(str, row['prefix'])))}</td>"
        f"<td>{row['input_count']}</td><td>{row['context_count']}</td></tr>"
        for row in audit["available_level_counts"]
    )
    count_table = (
        (
            "<h2>Available LIN depths</h2><table><thead><tr><th>Level</th><th>ST</th>"
            "<th>Prefix</th><th>Input</th><th>Context</th></tr></thead><tbody>"
            + counts
            + "</tbody></table>"
        )
        if counts
        else ""
    )
    report_path = output / "index.html"
    report_path.write_text(
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>ChronoClade · cgMLST</title><style>{report_styles()}</style></head>"
        '<body><main class="shell"><header class="identity"><div class="identity-mark">CHRONOCLADE · cgMLST</div><div class="identity-copy">'
        '<h1>Profile analysis blocks</h1></div></header><section class="stage"><div class="stage-body"><h2>Prepared dataset</h2>'
        f"<p>{records_count} samples; {len(products)} disjoint blocks. Status: {status}.</p>"
        "<p>Context counts describe compatible profiles in the prepared dataset, before assembly subsampling. "
        "Public retrieval provenance is saved with that dataset. Context assemblies are acquired only by the tree stage.</p>"
        f"<p>LIN depth: {lin_level}; HierCC: {html.escape(hiercc_level or 'Not selected')}.</p>"
        "<table><thead><tr><th>Block</th><th>Input</th><th>Context</th><th>Outputs</th></tr></thead>"
        f"<tbody>{items}</tbody></table>{count_table}<h2>Unresolved inputs</h2><ul>{unresolved or '<li>None</li>'}</ul>"
        '<a href="partitions.json">Partition counts and exclusions</a></div></section></main></body></html>'
    )
    return report_path
