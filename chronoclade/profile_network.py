"""Undirected exploratory location changes on a fixed profile tree.

Root agreement describes sensitivity to rooting, never transmission probability.
Nearest-country comparisons are observations and remain separate from reconstruction.
"""
from __future__ import annotations

import copy
import math
from collections import Counter

from chronoclade.context_refinement import _compatible, _lineage

_INPUT = {"local", "query", "focal", "input"}


def _country(row):
    value = str(row.get("country") or "").strip()
    return None if value.casefold() in {"", "unknown", "none", "na", "n/a"} else value


def _roots(tree, rows, limit=5):
    """Seed with an input; prioritise new countries, then farthest tree tips."""
    by_id = {r["sample_id"]: r for r in rows}
    candidates = sorted(t.name for t in tree.get_terminals() if t.name in by_id)
    if not candidates:
        return []
    inputs = [name for name in candidates if by_id[name].get("origin") in _INPUT]
    selected = [inputs[0] if inputs else candidates[0]]
    while len(selected) < min(limit, len(candidates)):
        seen = {_country(by_id[name]) for name in selected} - {None}
        remaining = [name for name in candidates if name not in selected]
        novel = [name for name in remaining if _country(by_id[name]) not in seen
                 and _country(by_id[name]) is not None]
        choices = novel or remaining
        # Stable tie-breaking by accession; actual selection is tree-distance based.
        selected.append(max(choices, key=lambda name: (
            min(tree.distance(name, root) for root in selected), -candidates.index(name))))
    return selected


def _possible_pairs(tree, countries, locations, root):
    rooted = copy.deepcopy(tree)
    if len(rooted.get_terminals()) > 1:
        rooted.root_with_outgroup(root)
    costs = {}
    for node in rooted.find_clades(order="postorder"):
        if node.is_terminal():
            observed = countries[node.name]
            costs[node] = {state: 0 if observed in {state, None} else math.inf
                           for state in locations}
        else:
            costs[node] = {
                state: sum(min(costs[child][target] + (state != target)
                               for target in locations) for child in node.clades)
                for state in locations}
    optimum = min(costs[rooted.root].values())
    states = {rooted.root: {s for s in locations if costs[rooted.root][s] == optimum}}
    possible, ambiguous_pairs, ambiguous = set(), set(), 0
    for parent in rooted.find_clades(order="preorder"):
        for child in parent.clades:
            allowed = set()
            for source in states[parent]:
                best = min(costs[child][target] + (source != target) for target in locations)
                allowed.update((source, target) for target in locations
                               if costs[child][target] + (source != target) == best)
            states[child] = {target for _, target in allowed}
            ambiguous += len(allowed) > 1
            pairs = {tuple(sorted((a, b))) for a, b in allowed if a != b}
            possible.update(pairs)
            if len(allowed) > 1:
                ambiguous_pairs.update(pairs)
    return possible, ambiguous_pairs, ambiguous


def _nearest(rows, comparisons):
    by_id = {row["sample_id"]: row for row in rows}
    groups = {}
    seen = set()
    for item in comparisons or []:
        qid, cid = item.get("query_id"), item.get("context_id")
        if qid not in by_id or not cid or (qid, cid) in seen:
            continue
        if item.get("status") != "matched":
            continue
        source = _country(by_id[qid]) or "Unknown"
        target = _country(by_id[cid] if cid in by_id else item) or "Unknown"
        seen.add((qid, cid))
        key = (source, target)  # Query country first, including same-country matches.
        groups.setdefault(key, []).append(item)
    result = []
    for (source, target), items in sorted(groups.items()):
        qs = sorted({i["query_id"] for i in items})
        cs = sorted({i["context_id"] for i in items})
        row = dict(source=source, target=target, query_ids=qs, context_ids=cs,
                   query_count=len(qs), context_count=len(cs), comparison_count=len(items))
        for field, alias in (("allele_differences", "allele_mismatches"),
                             ("shared_called_loci", "shared_loci")):
            values = [i[field] for i in items if i.get(field) is not None]
            row[alias + "_min"] = min(values) if values else None
            row[alias + "_max"] = max(values) if values else None
        result.append(row)
    return result


def _network(tree, rows, nearest):
    tip_ids = {tip.name for tip in tree.get_terminals()}
    rows = [r for r in rows if r["sample_id"] in tip_ids]
    countries = {r["sample_id"]: _country(r) for r in rows}
    locations = sorted(set(countries.values()) - {None})
    roots = _roots(tree, rows) if locations else []
    counts, ambiguity_counts, ambiguity = Counter(), Counter(), 0
    for root in roots:
        possible, ambiguous_pairs, ambiguous = _possible_pairs(tree, countries, locations, root)
        counts.update(possible)
        ambiguity_counts.update(ambiguous_pairs)
        ambiguity += ambiguous
    nodes = []
    for country in locations:
        subset = [r for r in rows if _country(r) == country]
        inputs = sum(r.get("origin") in _INPUT for r in subset)
        nodes.append(dict(country=country, count=len(subset), input_count=inputs,
                          context_count=len(subset) - inputs, is_input_country=bool(inputs)))
    return dict(
        method="undirected possible maximum-parsimony changes across diverse tip roots",
        interpretation="Possible location changes, not demonstrated transmission. Root fractions "
        "measure sensitivity on a fixed tree, not probability or confidence; sampling and tree "
        "uncertainty are not estimated.",
        root_selection="input seed, then country diversity and farthest tree tips",
        uncertainty_scope="fixed tree; no topology or sampling uncertainty assessment",
        tested_roots=roots, root_count=len(roots),
        ambiguous_edge_reconstructions=ambiguity, nodes=nodes,
        unknown_country_count=sum(value is None for value in countries.values()),
        input_countries=[n["country"] for n in nodes if n["is_input_country"]],
        input_ids=sorted(r["sample_id"] for r in rows if r.get("origin") in _INPUT),
        sample_count=len(rows),
        edges=[dict(source=a, target=b, roots_with_possible_change=count,
                    roots_tested=len(roots), root_fraction=count / len(roots), uncertain=True,
                    roots_with_ambiguous_change=ambiguity_counts[(a, b)])
               for (a, b), count in sorted(counts.items())],
        nearest_edges=_nearest(rows, nearest))


def build_profile_network(tree, rows, nearest_neighbours=None, subgroup_depth=5):
    """Build full-pool network and exact, version-compatible input LIN-prefix views."""
    network = _network(tree, rows, nearest_neighbours)
    groups = {}
    for row in rows:
        lineage = _lineage(row, "cglin", subgroup_depth)
        if row.get("origin") in _INPUT and lineage is not None:
            groups.setdefault(lineage, []).append(row)
    views = []
    for index, (lineage, inputs) in enumerate(sorted(groups.items()), 1):
        subset = [row for row in rows if _lineage(row, "cglin", subgroup_depth) == lineage
                  and _compatible(inputs[0], row, "cglin")]
        ids = {r["sample_id"] for r in subset}
        pruned = copy.deepcopy(tree)
        for tip in list(pruned.get_terminals()):
            if tip.name not in ids:
                pruned.prune(tip)
        data = _network(pruned, subset, nearest_neighbours)
        input_ids = set(data["input_ids"])
        subgroup_nearest = [item for item in nearest_neighbours or []
                            if item.get("query_id") in input_ids]
        data["nearest_edges"] = _nearest(rows, subgroup_nearest)
        data["nearest_search_scope"] = "complete selected CG pool; subgroup input queries"
        prefix = ",".join(map(str, lineage[-1]))
        views.append(dict(id=f"lin_{index}", label=f"LIN level {subgroup_depth}: {prefix} ({len(inputs)} inputs)",
                          scope=dict(kind="cglin", level=subgroup_depth, prefix=prefix,
                                     scheme=lineage[0], version=lineage[1]), **data))
    network["views"] = views
    return network


def draw_profile_network(network, path, focus_inputs=True):
    """Export an undirected figure, emphasising links touching input countries."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    inputs = set(network["input_countries"])
    edges = [e for e in network["edges"] if not focus_inputs or
             e["source"] in inputs or e["target"] in inputs]
    visible = inputs | {e[k] for e in edges for k in ("source", "target")}
    nodes = [n for n in network["nodes"] if not focus_inputs or n["country"] in visible]
    # Input countries occupy the centre; comparisons around them keep focus plots legible.
    public = [n for n in nodes if n["country"] not in inputs]
    focal = [n for n in nodes if n["country"] in inputs]
    positions = {}
    for group, radius in ((public, 1.0), (focal, 0.25 if public else 0.7)):
        for index, node in enumerate(group):
            angle = 2 * math.pi * index / max(1, len(group))
            positions[node["country"]] = (radius * math.cos(angle), radius * math.sin(angle))
    fig, ax = plt.subplots(figsize=(10, 8))
    for edge in edges:
        a, b = positions[edge["source"]], positions[edge["target"]]
        fraction = edge["root_fraction"]
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#777777", alpha=0.55,
                linewidth=0.7 + fraction, linestyle="-" if fraction == 1 and not edge.get("roots_with_ambiguous_change")
                else "--", zorder=1)
    for node in nodes:
        x, y = positions[node["country"]]
        ax.scatter(x, y, s=45 + 35 * math.sqrt(node["count"]),
                   color="#2166ac" if node["is_input_country"] else "#bbbbbb", zorder=2)
        ax.annotate(f"{node['country']} ({node['count']})", (x, y), xytext=(5, 6),
                    textcoords="offset points", fontsize=8, zorder=3)
    if not nodes:
        ax.text(0.5, 0.5, "No location links available for input countries", transform=ax.transAxes,
                ha="center")
    ax.legend(handles=[Line2D([], [], marker="o", color="#2166ac", linestyle="",
                             label="Country containing input genomes"),
                       Line2D([], [], color="#777777", label="Possible across roots, no detected state ambiguity"),
                       Line2D([], [], color="#777777", linestyle="--",
                              label="Root-sensitive or ambiguous location assignment")],
              loc="upper center", bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=8)
    ax.set_title("Possible country connections · " +
                 ("input-country focus" if focus_inputs else "all countries"))
    ax.text(0.5, -0.13, "Node size: sampled genomes. Line style: root sensitivity, not confidence.\n"
            "Undirected links are not demonstrated transmission.", transform=ax.transAxes,
            ha="center", fontsize=8)
    ax.set(xlim=(-1.3, 1.5), ylim=(-1.3, 1.3))
    ax.axis("off")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return str(path)
