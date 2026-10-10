"""Independent ESM2 exploration report using shared offline report components."""
from html import escape
from pathlib import Path
from chronoclade.report_components.styles import report_styles


def write_genome_report(summary, directory, *, ordination_html, ensemble_path=None,
                        baseline_report=None, temporal_report=None):
    directory = Path(directory)
    counts, panel = summary["counts"], summary["panel"]
    def link(path, title):
        return f'<p><a href="{escape(Path(path).relative_to(directory).as_posix())}">{escape(title)}</a></p>' if path else ""
    rows = "".join('<tr>' + ''.join(f'<td>{escape(str(value))}</td>' for value in
                   (sample["sample_id"], sample.get("role"), sample["status"], sample.get("reason") or "Included",
                    ", ".join(sample.get("missing_panel_loci", [])))) + '</tr>' for sample in summary["samples"])
    neighbours = "".join(f'<tr><td>{escape(row["sample_id"])}</td><td>' + ', '.join(
        escape(item["sample_id"]) + f' ({item["distance"]:.5g})' for item in row["neighbours"]) + '</td></tr>' for row in summary["neighbours"])
    body = f'''<header class="identity"><div class="identity-mark">CHRONOCLADE · ESM2</div>
<div class="identity-copy"><h1>Genome protein exploration</h1><p>Experimental locus-preserving embeddings on frozen prepared samples.</p></div></header>
<nav class="contents"><a href="#scope">Scope</a><a href="#ordination">Ordination</a><a href="#neighbours">Neighbours</a><a href="#coverage">Coverage</a><a href="#evidence">Evidence</a></nav>
<section class="stage" id="scope"><div class="stage-body"><h2>Fixed protein panel</h2>
<p>{counts["total_samples"]} prepared genomes; {counts["eligible_samples"]} eligible genomes; {panel["count"]} common loci; {len(panel["variable_protein_loci"])} variable protein loci.</p>
<p>Distance is the equal-weight mean of locus-specific cosine distances. Each locus retains its own protein identity. Missing vectors are excluded and never replaced by zeros.</p>
<p class="caveat">This route remains experimental. Synonymous DNA changes can produce identical proteins. Embedding groups are exploratory connected components and have no calibrated biological threshold. No embedding ancestor, substitution rate, clock or TMRCA is inferred.</p>
{link(baseline_report, "Conventional cgMLST guide tree and biological location network")}
{link(temporal_report, "Descriptive collection-date comparisons and residuals")}</div></section>
<section class="stage" id="ordination"><div class="stage-body"><h2>Genome ordination</h2>
<p>Coordinates use protein chord distance: sqrt(2 × mean locus cosine distance). This preserves locus identity and locates every eligible genome. Country, collection date, role and group controls use the shared report viewer.</p>{ordination_html}</div></section>
<section class="stage" id="neighbours"><div class="stage-body"><h2>Nearest genomes</h2><p>All ties at the requested boundary are retained. Values are mean locus cosine distances.</p><div class="table-scroll"><table><thead><tr><th>Genome</th><th>Neighbours (embedding distance)</th></tr></thead><tbody>{neighbours}</tbody></table></div></div></section>
<section class="stage" id="coverage"><div class="stage-body"><h2>Mapping and panel coverage</h2><p>Panel selection: {escape(panel["selection"])}. {len(panel["omitted_scheme_loci"])} scheme loci are omitted; exact IDs and mapped masks are downloadable.</p><div class="table-scroll"><table><thead><tr><th>Genome</th><th>Role</th><th>Status</th><th>Reason</th><th>Missing panel loci</th></tr></thead><tbody>{rows}</tbody></table></div></div></section>
<section class="stage" id="evidence"><div class="stage-body"><h2>Evidence and selections</h2>
<p>Allele neighbour agreement uses the same genomes and frozen loci. Corrected SNP validation has not been supplied; runtime success does not validate within-lineage biological performance.</p>
{link(ensemble_path, "Validated selection ensemble for the independent tree stage")}
<ul class="downloads"><li><a href="genome_analysis.json" download>Analysis, exclusions, model and locus evidence (JSON)</a></li><li><a href="genome_vectors.npz" download>Distances, ordination, locus protein IDs and missing masks (NPZ)</a></li><li><a href="pairwise_distances.csv" download>Genome distances and allele/protein comparison (CSV)</a></li></ul></div></section>'''
    html = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ChronoClade · Genome protein exploration</title><style>' + report_styles() + '</style></head><body><main class="shell">' + body + '</main></body></html>'
    path = directory / "index.html"
    path.write_text(html, encoding="utf-8")
    return path


def write_genome_index(products, audit, output):
    items = "".join(f'<li><a href="{escape(row["report"])}">{escape(row["block_id"])}</a> · {row["eligible_genomes"]} eligible genomes · {escape(row["status"])}</li>' for row in products)
    detail = escape(str(audit.get("unresolved_inputs", [])))
    path = Path(output) / "index.html"
    path.write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ChronoClade · ESM2 blocks</title><style>' + report_styles() + '</style></head><body><main class="shell"><header class="identity"><h1>ESM2 genome exploration</h1></header><section class="stage"><div class="stage-body"><h2>Compatible frozen lineage blocks</h2><p>Same explicit lineage partition policy as cgMLST; available frozen context only. Embeddings remain experimental.</p><ul>' + items + '</ul><p>Unresolved inputs: ' + detail + '</p></div></section></main></body></html>', encoding="utf-8")
    return path
