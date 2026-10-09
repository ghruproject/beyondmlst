"""Closest analysed context genomes and genetic relationships, independent of dating."""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import shutil
from collections import Counter
from html import escape
from pathlib import Path

from Bio import Phylo

from chronoclade.errors import WorkflowError
from chronoclade.metadata import Sample


class NeighbourhoodError(WorkflowError):
    """Neighbourhood inputs cannot establish trustworthy sample identities."""


_SCOPE = (
    "Only genomes in this analysis were compared. These are the closest relatives within "
    "the analysed context set, not necessarily the closest in all public data. "
    "Similarity alone does not establish transmission, acquisition country or migration direction."
)
_FIELDS = [
    "focal_sample",
    "focal_country",
    "focal_location",
    "focal_date",
    "context_sample",
    "context_country",
    "context_location",
    "context_date",
    "context_source",
    "source_genome_id",
    "public_provenance",
    "metric",
    "rank",
    "distance",
    "clonal_snps",
    "callable_sites",
    "snp_proportion",
    "patristic_distance",
    "tied_at_distance",
    "scope",
]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _context(sample: Sample) -> bool:
    return sample.origin.strip().casefold() == "context"


def _source_rows(directory: Path) -> dict[str, dict]:
    manifest = directory / "context_manifest.tsv"
    if not manifest.is_file():
        return {}
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    for row in rows:
        ident = row.get("sample_id", "")
        if ident in result and result[ident] != row:
            raise NeighbourhoodError("Conflicting context manifest sample IDs")
        result[ident] = row
    return result


def _pairs(path: Path, sample_ids: set[str]) -> tuple[dict, list[dict]]:
    pairs, unavailable = {}, []
    if not path.is_file():
        return pairs, [{"reason": "SNP table unavailable"}]
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_1", "sample_2", "clonal_snps", "callable_sites"}
        if not required.issubset(reader.fieldnames or []):
            return pairs, [{"reason": "SNP table is missing required columns"}]
        for index, row in enumerate(reader, 2):
            ids = (row["sample_1"], row["sample_2"])
            if ids[0] not in sample_ids or ids[1] not in sample_ids:
                unavailable.append({"row": index, "reason": "pair contains an unanalysed sample"})
                continue
            if ids[0] == ids[1]:
                unavailable.append({"row": index, "reason": "self comparison excluded"})
                continue
            key = tuple(sorted(ids))
            if key in pairs:
                raise NeighbourhoodError("SNP table contains duplicate sample-pair records")
            try:
                snps = int(row["clonal_snps"])
                callable_sites = int(row["callable_sites"])
                if snps < 0 or callable_sites <= 0 or snps > callable_sites:
                    raise ValueError
            except (TypeError, ValueError):
                unavailable.append(
                    {
                        "row": index,
                        "sample_1": ids[0],
                        "sample_2": ids[1],
                        "reason": "invalid SNP count or missing/zero callable sites",
                    }
                )
                continue
            pairs[key] = {
                "clonal_snps": snps,
                "callable_sites": callable_sites,
                "snp_proportion": snps / callable_sites,
            }
    return pairs, unavailable


def _read_tree(path: Path):
    try:
        tree = Phylo.read(path, "newick")
    except (OSError, ValueError, IndexError):
        return None, "Genetic tree is missing or unreadable"
    ids = [tip.name for tip in tree.get_terminals()]
    if any(not ident for ident in ids) or len(ids) != len(set(ids)):
        return None, "Genetic tree has missing or duplicate tip IDs"
    for clade in tree.find_clades():
        value = clade.branch_length
        if value is None and clade is tree.root:
            continue  # The edge above the root is not a between-tip genetic distance.
        if value is None or not math.isfinite(value) or value < 0:
            return None, "Genetic tree has missing, negative or non-finite branch lengths"
    return tree, ""


def _ranking(comparisons: list[dict], metric: str, top_n: int) -> list[dict]:
    field = "clonal_snps" if metric == "clonal_snps" else "patristic_distance"
    ordered = sorted(comparisons, key=lambda row: (row[field], row["context_sample"]))
    if not ordered:
        return []
    cutoff = ordered[min(top_n, len(ordered)) - 1][field]
    frequencies = Counter(row[field] for row in ordered)
    first_rank = {}
    for index, row in enumerate(ordered, 1):
        first_rank.setdefault(row[field], index)
    return [
        {
            **row,
            "metric": metric,
            "rank": first_rank[row[field]],
            "distance": row[field],
            "tied_at_distance": frequencies[row[field]],
            "scope": _SCOPE,
        }
        for row in ordered
        if row[field] <= cutoff
    ]


def _render_page(
    tree,
    sample_by_id: dict[str, Sample],
    page_ids: list[str],
    path: Path,
    *,
    page: int,
    page_count: int,
    total: int,
) -> None:
    """Draw real branch distances; avoid Bio.Phylo's unit-length fallback for zero trees."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    selected = set(page_ids)
    view = copy.deepcopy(tree)
    for tip in list(view.get_terminals()):
        if tip.name not in selected:
            view.prune(tip)
    tips = view.get_terminals()
    positions = {tip: float(len(tips) - index - 1) for index, tip in enumerate(tips)}
    depths = {view.root: 0.0}

    def locate(clade):
        if clade.is_terminal():
            return positions[clade]
        for child in clade.clades:
            depths[child] = depths[clade] + float(child.branch_length or 0.0)
            locate(child)
        positions[clade] = (
            min(positions[c] for c in clade.clades) + max(positions[c] for c in clade.clades)
        ) / 2
        return positions[clade]

    locate(view.root)
    maximum = max(depths.values(), default=0.0)
    label_pad = maximum * 0.025 if maximum else 0.0025
    height = max(4.8, len(tips) * 0.62 + 1.4)
    with matplotlib.rc_context({"svg.fonttype": "none", "text.usetex": False}):
        figure, axis = plt.subplots(figsize=(8, height))
        figure.subplots_adjust(left=0.07, right=0.43, top=1 - 0.8 / height, bottom=0.65 / height)
        for clade in view.find_clades(order="preorder"):
            if clade.clades:
                values = [positions[c] for c in clade.clades]
                axis.plot([depths[clade]] * 2, [min(values), max(values)], color="#17191f", lw=0.8)
            for child in clade.clades:
                axis.plot(
                    [depths[clade], depths[child]], [positions[child]] * 2, color="#17191f", lw=0.8
                )
        for tip in tips:
            sample = sample_by_id[tip.name]
            focal = not _context(sample)
            location = sample.location or "Unknown location"
            day = sample.collection_date or "Unknown date"
            role = "FOCAL" if focal else "context"
            # Separate ID from metadata so the inline view keeps readable type.
            label = f"{role} · {tip.name}\n{location} · {day}"
            axis.scatter(
                depths[tip],
                positions[tip],
                s=28 if focal else 16,
                color="#2855a6" if focal else "#5f6470",
                zorder=3,
            )
            axis.text(
                depths[tip] + label_pad,
                positions[tip],
                label,
                va="center",
                fontsize=14,
                color="#2855a6" if focal else "#5f6470",
                fontweight="bold" if focal else "normal",
                clip_on=False,
                parse_math=False,
                linespacing=1.2,
            )
        axis.set_ylim(-0.6, max(1, len(tips) - 0.4))
        axis.set_xlim(0, maximum + label_pad * 2 if maximum else 0.1)
        axis.set_yticks([])
        axis.set_xlabel("Substitutions per site", fontsize=12, labelpad=8)
        axis.tick_params(axis="x", labelsize=12)
        axis.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3), useMathText=False)
        for side in ("left", "top", "right"):
            axis.spines[side].set_visible(False)
        figure.suptitle(
            "Genetic relationships before dating",
            x=0.07,
            y=1 - 0.08 / height,
            ha="left",
            fontweight="bold",
            fontsize=15,
        )
        note = (
            f"Page {page}/{page_count} · {len(tips)} of {total} tips · focal samples bold and blue"
            if page_count > 1
            else f"All {total} analysed tips · focal samples bold and blue"
        )
        figure.text(0.07, 1 - 0.43 / height, note, fontsize=11, ha="left")
        if maximum == 0:
            axis.text(
                0.01,
                0.02,
                "No positive branch span in this view",
                transform=axis.transAxes,
                fontsize=11,
            )
        figure.savefig(
            path, format="svg", bbox_inches="tight", facecolor="white", metadata={"Date": None}
        )
        figure.savefig(
            path.with_suffix(".png"),
            dpi=160,
            bbox_inches="tight",
            facecolor="white",
            metadata={"Software": "ChronoClade"},
        )
        plt.close(figure)


def build_neighbourhood_evidence(
    *, tree: Path, samples: list[Sample], output: Path, top_n: int = 3, max_tips_per_page: int = 36
) -> dict:
    """Rank final clonal SNP and tree distances separately, preserving every cutoff tie."""
    if top_n < 1 or max_tips_per_page < 2:
        raise NeighbourhoodError(
            "Ranking count must be positive and tree pages need at least two tips"
        )
    sample_by_id = {sample.sample_id: sample for sample in samples}
    if len(sample_by_id) != len(samples) or not samples:
        raise NeighbourhoodError("Neighbourhood requires unique nonempty sample IDs")
    output.mkdir(parents=True, exist_ok=True)
    pairwise = output / "clonal_pairwise_distances.tsv"
    pairs, invalid_pairs = _pairs(pairwise, set(sample_by_id))
    sources = _source_rows(output)
    contexts = sorted((s for s in samples if _context(s)), key=lambda s: s.sample_id)
    focal = sorted((s for s in samples if not _context(s)), key=lambda s: s.sample_id)
    phylogeny, tree_error = _read_tree(tree)
    tree_ids = {tip.name for tip in phylogeny.get_terminals()} if phylogeny else set()
    missing = sorted(set(sample_by_id) - tree_ids) if phylogeny else sorted(sample_by_id)
    extra = sorted(tree_ids - set(sample_by_id))
    rows, summaries = [], []
    for sample in focal:
        snp_comparisons, tree_comparisons = [], []
        for context in contexts:
            source = sources.get(context.sample_id, {})
            provenance = (
                "verified_provider"
                if source.get("source") == "pathogenwatch"
                else "unverified_supplied_context"
            )
            record = {
                "focal_sample": sample.sample_id,
                "focal_country": "",
                "focal_location": sample.location,
                "focal_date": sample.collection_date,
                "context_sample": context.sample_id,
                "context_country": source.get("country", ""),
                "context_location": context.location or "Unknown",
                "context_date": context.collection_date,
                "context_source": source.get("source") or "supplied_context",
                "source_genome_id": source.get("source_genome_id")
                or source.get("sample_accession")
                or "",
                "public_provenance": provenance,
                "clonal_snps": None,
                "callable_sites": None,
                "snp_proportion": None,
                "patristic_distance": None,
            }
            pair = pairs.get(tuple(sorted((sample.sample_id, context.sample_id))))
            if pair:
                record.update(pair)
            if phylogeny and sample.sample_id in tree_ids and context.sample_id in tree_ids:
                distance = float(phylogeny.distance(sample.sample_id, context.sample_id))
                record["patristic_distance"] = distance
                if math.isfinite(distance) and distance >= 0:
                    tree_comparisons.append(dict(record))
                else:
                    record["patristic_distance"] = None
            if pair:
                snp_comparisons.append(dict(record))
        ranked_snp = _ranking(snp_comparisons, "clonal_snps", top_n)
        ranked_tree = _ranking(tree_comparisons, "patristic_distance", top_n)
        rows.extend(ranked_snp)
        rows.extend(ranked_tree)
        best_snp = [r["context_sample"] for r in ranked_snp if r["rank"] == 1]
        best_tree = [r["context_sample"] for r in ranked_tree if r["rank"] == 1]
        summaries.append(
            {
                "focal_sample": sample.sample_id,
                "status": "no_context"
                if not contexts
                else "available"
                if ranked_snp or ranked_tree
                else "unavailable",
                "snp_status": "available" if ranked_snp else "unavailable",
                "tree_status": "available" if ranked_tree else "unavailable",
                "snp_comparisons": len(snp_comparisons),
                "tree_comparisons": len(tree_comparisons),
                "nearest_by_clonal_snps": best_snp,
                "nearest_by_tree_distance": best_tree,
                "rankings_agree": set(best_snp) == set(best_tree)
                if best_snp and best_tree
                else None,
                "missing_snp_context_ids": [
                    c.sample_id
                    for c in contexts
                    if tuple(sorted((sample.sample_id, c.sample_id))) not in pairs
                ],
                "missing_tree_context_ids": [
                    c.sample_id for c in contexts if c.sample_id not in tree_ids
                ],
                "focal_missing_from_tree": sample.sample_id not in tree_ids,
            }
        )
    pages = []
    for alias in ("genetic_tree.svg", "genetic_tree.png"):
        (output / alias).unlink(missing_ok=True)
    # Remove artifacts from a previous run before writing current figures.
    for previous in output.glob("genetic_relationship_tree*.svg"):
        previous.unlink()
    for previous in output.glob("genetic_relationship_tree*.png"):
        previous.unlink()
    if phylogeny:
        ids = [tip.name for tip in phylogeny.get_terminals() if tip.name in sample_by_id]
        count = math.ceil(len(ids) / max_tips_per_page)
        for index in range(count):
            page_ids = ids[index * max_tips_per_page : (index + 1) * max_tips_per_page]
            name = (
                "genetic_relationship_tree.svg"
                if index == 0
                else f"genetic_relationship_tree_page{index + 1}.svg"
            )
            _render_page(
                phylogeny,
                sample_by_id,
                page_ids,
                output / name,
                page=index + 1,
                page_count=count,
                total=len(ids),
            )
            pages.append(
                {
                    "page": index + 1,
                    "svg": name,
                    "png": name.replace(".svg", ".png"),
                    "tip_ids": page_ids,
                    "omitted_tip_count": len(ids) - len(page_ids),
                }
            )
        (output / "genetic_relationship_tree.newick").write_bytes(tree.read_bytes())
        if pages:
            shutil.copyfile(output / pages[0]["svg"], output / "genetic_tree.svg")
            shutil.copyfile(output / pages[0]["png"], output / "genetic_tree.png")
    else:
        (output / "genetic_relationship_tree.newick").unlink(missing_ok=True)
    result = {
        "schema_version": 1,
        "scope": _SCOPE,
        "top_n_before_cutoff_ties": top_n,
        "context_genomes": len(contexts),
        "focal_genomes": len(focal),
        "rows": rows,
        "focal_summaries": summaries,
        "invalid_snp_pairs": invalid_pairs,
        "tree_status": "available" if phylogeny else "unavailable",
        "tree_message": tree_error,
        "missing_tree_tip_ids": missing,
        "extra_tree_tip_ids": extra,
        "tree_pages": pages,
        "source_tree_sha256": _sha(tree) if tree.is_file() else None,
        "source_snp_table_sha256": _sha(pairwise) if pairwise.is_file() else None,
        "source_context_manifest_sha256": _sha(output / "context_manifest.tsv")
        if (output / "context_manifest.tsv").is_file()
        else None,
        "sample_metadata_sha256": hashlib.sha256(
            json.dumps(
                [
                    (s.sample_id, s.location, s.collection_date, s.origin)
                    for s in sorted(samples, key=lambda s: s.sample_id)
                ],
                separators=(",", ":"),
            ).encode()
        ).hexdigest(),
    }
    result["outputs"] = {
        "nearest_neighbours": str(output / "nearest_neighbours.tsv"),
        "nearest_neighbours_json": str(output / "nearest_neighbours.json"),
        "genetic_tree": str(output / "genetic_tree.svg") if pages else None,
        "genetic_tree_png": str(output / "genetic_tree.png") if pages else None,
    }
    for suffix, delimiter in (("tsv", "\t"), ("csv", ",")):
        with (output / f"nearest_neighbours.{suffix}").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=_FIELDS, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
    (output / "nearest_neighbours.json").write_text(
        json.dumps(result, sort_keys=True, indent=2), encoding="utf-8"
    )
    (output / "genetic_relationship_tree.json").write_text(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "tree_status",
                    "tree_message",
                    "missing_tree_tip_ids",
                    "extra_tree_tip_ids",
                    "tree_pages",
                    "source_tree_sha256",
                )
            },
            sort_keys=True,
            indent=2,
        ),
        encoding="utf-8",
    )
    return result


def neighbourhood_report_html(directory: Path) -> str:
    """Tree first, nearest-SNP summaries next, with full rankings available on demand."""
    path = directory / "nearest_neighbours.json"
    if not path.is_file():
        return ""
    data = json.loads(path.read_text(encoding="utf-8"))
    figures = []
    for page in data["tree_pages"]:
        suffix = f" (page {page['page']})" if len(data["tree_pages"]) > 1 else ""
        figures.append(
            f'<figure><a href="{escape(page["svg"], quote=True)}"><img src="{escape(page["svg"], quote=True)}" '
            f'alt="Genetic relationships before dating{suffix}; focal samples emphasised"></a>'
            f"<figcaption>Genetic relationships{suffix}: {len(page['tip_ids'])} tips shown; "
            f"{page['omitted_tip_count']} omitted from this view. Branches measure substitutions per site, "
            "not time. Labels give recorded locations and collection dates; ancestral locations are not inferred. "
            f'<a href="{escape(page["svg"], quote=True)}">Open full-resolution tree</a> · '
            f'<a href="{escape(page["png"], quote=True)}">Download PNG</a>.</figcaption></figure>'
        )
    tree_visual = (
        figures[0]
        if figures
        else f"<p>{escape(data['tree_message'] or 'Genetic tree unavailable.')}</p>"
    )
    if len(figures) > 1:
        tree_visual += (
            "<details><summary>All genetic tree pages</summary>"
            + "".join(figures[1:])
            + "</details>"
        )
    if data["missing_tree_tip_ids"]:
        tree_visual += (
            f"<p>{len(data['missing_tree_tip_ids'])} analysed samples are missing from the tree; "
            "their SNP comparisons remain available. Missing IDs are in the JSON audit.</p>"
        )
    summaries, ranked_tables = [], []
    for summary in data["focal_summaries"]:
        focal_id = summary["focal_sample"]
        selected = [row for row in data["rows"] if row["focal_sample"] == focal_id]
        closest = [row for row in selected if row["metric"] == "clonal_snps" and row["rank"] == 1]
        agreement = (
            "The closest sets differ between SNP and tree distance; review both measures."
            if summary["rankings_agree"] is False
            else "The genetic tree identifies the same closest set."
            if summary["rankings_agree"] is True
            else "A comparison measure is unavailable; agreement cannot be assessed."
        )
        title = f"<h3>Closest relatives of {escape(focal_id)}</h3>"
        if closest:
            body = "".join(
                "<tr>"
                + "".join(
                    f'<td data-label="{escape(label, quote=True)}">{escape(str(value))}</td>'
                    for label, value in zip(
                        (
                            "Context genome",
                            "Recorded location",
                            "Collection date",
                            "DNA differences",
                            "Comparable sites",
                        ),
                        (
                            row["context_sample"],
                            row["context_country"] or row["context_location"],
                            row["context_date"] or "Unknown",
                            row["clonal_snps"],
                            f"{row['callable_sites']:,}",
                        ),
                    )
                )
                + "</tr>"
                for row in closest
            )
            count = summary["snp_comparisons"]
            summaries.append(
                title
                + f"<p>Compared {count} of {data['context_genomes']} analysed context genomes; "
                f"{len(closest)} closest {'genome' if len(closest) == 1 else 'genomes'} by final clonal SNPs. "
                f"{agreement}</p>" + '<div class="table-scroll" tabindex="0" role="region" '
                'aria-label="Closest relatives by final clonal SNPs"><table><thead><tr>'
                "<th>Context genome</th><th>Recorded location</th><th>Date</th><th>DNA differences (SNPs)</th>"
                "<th>Comparable sites</th></tr></thead><tbody>" + body + "</tbody></table></div>"
            )
        else:
            message = (
                "No context genomes were analysed."
                if summary["status"] == "no_context"
                else "Final SNP comparisons are unavailable; review callable-site coverage in the audit."
            )
            summaries.append(title + "<p>" + message + " " + agreement + "</p>")
        if selected:
            body = []
            for row in selected:
                values = (
                    "Final clonal SNPs"
                    if row["metric"] == "clonal_snps"
                    else "Genetic tree distance",
                    row["rank"],
                    row["context_sample"],
                    row["context_country"] or row["context_location"],
                    row["context_date"] or "Unknown",
                    "—" if row["clonal_snps"] is None else row["clonal_snps"],
                    "—" if row["callable_sites"] is None else f"{row['callable_sites']:,}",
                    "—" if row["snp_proportion"] is None else f"{row['snp_proportion']:.4g}",
                    "—"
                    if row["patristic_distance"] is None
                    else f"{row['patristic_distance']:.4g}",
                    row["tied_at_distance"],
                )
                body.append(
                    "<tr>" + "".join(f"<td>{escape(str(v))}</td>" for v in values) + "</tr>"
                )
            ranked_tables.append(
                f'<h4>{escape(focal_id)}</h4><div class="table-scroll" tabindex="0" '
                'role="region" aria-label="All ranked context relatives"><table><thead><tr>'
                "<th>Ranked by</th><th>Rank</th><th>Context genome</th><th>Recorded location</th>"
                "<th>Date</th><th>Clonal SNPs</th><th>Callable sites</th><th>SNP/site</th>"
                "<th>Tree distance</th><th>Tied genomes</th></tr></thead><tbody>"
                + "".join(body)
                + "</tbody></table></div>"
            )
    details = (
        (
            "<details><summary>All ranked relatives and comparison measures</summary>"
            "<p>SKA distances screened candidate assemblies. Final clonal SNPs count differences in the "
            "recombination-filtered alignment. SNP/site uses each pair’s callable-site denominator. "
            "Tree distance sums branches between tips in substitutions per site, independently of dating. "
            "The rankings may differ; all exact ties at the top-rank cutoff are retained.</p>"
            + "".join(ranked_tables)
            + "</details>"
        )
        if ranked_tables
        else ""
    )
    return (
        '<div class="neighbourhood">'
        + tree_visual
        + "<p>"
        + escape(_SCOPE)
        + "</p><p>DNA differences (SNPs) are counted after removing inferred recombination. Comparable sites are the DNA positions that could be compared between the two genomes.</p>"
        + "".join(summaries)
        + details
        + '<p><a href="nearest_neighbours.tsv">Download nearest-relative TSV</a> · '
        '<a href="nearest_neighbours.csv">CSV</a> · <a href="nearest_neighbours.json">Audit JSON</a>'
        + (
            ' · <a href="genetic_relationship_tree.newick">Source genetic tree</a>'
            if data["tree_pages"]
            else ""
        )
        + "</p><p>Manually supplied context without a provider manifest has unverified public provenance. "
        "Genetic similarity alone does not establish transmission.</p></div>"
    )
