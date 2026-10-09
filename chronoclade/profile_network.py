"""Weighted undirected country changes and exact optimal-history count ranges.

Ranges condition on a fixed tree. Root agreement is a separate sensitivity audit,
never transmission probability.
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
            observed = countries.get(node.name)
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


def country_palette(countries):
    """Stable country colours independent of the selected view."""
    import colorsys
    import hashlib

    return {country: "#" + "".join(f"{round(v * 255):02x}" for v in
            colorsys.hsv_to_rgb(int.from_bytes(hashlib.sha256(country.encode()).digest()[:4],
                                             "big") / 2**32, 0.58, 0.72))
            for country in countries}


def _reconstruct(tree, countries, locations, seed=1729):
    """Primary Sankoff optimum and exact independent pair-count extrema.

    Secondary objectives only traverse primary-optimal child transitions. The
    seeded representative is a reproducible local tie choice, not a uniform
    sample of histories. Pair extrema need not be jointly attainable.
    """
    import random
    import numpy as np
    from io import StringIO
    from Bio import Phylo

    nodes = list(tree.find_clades(order="preorder"))
    ids = {node: i for i, node in enumerate(nodes)}
    parents = {child: node for node in nodes for child in node.clades}
    k = len(locations)
    pairs = [(a, b) for i, a in enumerate(locations) for b in locations[i + 1:]]
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    state_index = {state: i for i, state in enumerate(locations)}
    p = len(pairs)
    costs, lows, highs, allowed = {}, {}, {}, {}
    for node in reversed(nodes):
        costs[node] = np.zeros(k)
        if not k:
            continue
        if node.is_terminal():
            observed = countries.get(node.name)
            if observed is not None:
                costs[node][:] = np.inf
                costs[node][state_index[observed]] = 0
            continue
        lows[node] = np.zeros((k, p), dtype=np.int32)
        highs[node] = np.zeros((k, p), dtype=np.int32)
        for child in node.clades:
            child_cost = costs[child]
            best = np.minimum(child_cost, child_cost.min() + 1)
            costs[node] += best
            options = []
            for source in range(k):
                targets = np.flatnonzero(child_cost + (np.arange(k) != source) == best[source])
                options.append(targets)
                if child.is_terminal():
                    lower = np.zeros((len(targets), p), dtype=np.int32)
                    upper = lower.copy()
                else:
                    lower = lows[child][targets].copy()
                    upper = highs[child][targets].copy()
                for i, target in enumerate(targets):
                    if source != target:
                        col = pair_index[tuple(sorted((locations[source], locations[target])))]
                        lower[i, col] += 1
                        upper[i, col] += 1
                lows[node][source] += lower.min(axis=0)
                highs[node][source] += upper.max(axis=0)
            allowed[child] = options
        for child in node.clades:
            if not child.is_terminal():
                del lows[child], highs[child]
    optimum = int(costs[tree.root].min()) if k else 0
    root_states = np.flatnonzero(costs[tree.root] == optimum) if k else []
    minima = lows[tree.root][root_states].min(axis=0) if k and not tree.root.is_terminal() else []
    maxima = highs[tree.root][root_states].max(axis=0) if k and not tree.root.is_terminal() else []
    rng = random.Random(seed)
    assigned, global_states, counts = {}, {}, Counter()
    if k:
        assigned[tree.root] = int(rng.choice(list(root_states)))
        global_states[tree.root] = set(root_states)
        for node in nodes:
            for child in node.clades:
                source = assigned[node]
                target = int(rng.choice(list(allowed[child][source])))
                assigned[child] = target
                global_states[child] = {int(t) for s in global_states[node]
                                        for t in allowed[child][s]}
                if source != target:
                    counts[tuple(sorted((locations[source], locations[target])))] += 1
    audit_nodes, branches = [], []
    for node in nodes:
        state = locations[assigned[node]] if k else None
        observed = countries.get(node.name) if node.is_terminal() else None
        unknown = node.is_terminal() and observed is None
        audit_nodes.append(dict(id=ids[node], parent_id=ids[parents[node]] if node in parents else None,
                                children_ids=[ids[c] for c in node.clades], name=node.name,
                                is_tip=node.is_terminal(), state=None if unknown else state,
                                inferred_state=state, observed_country=observed,
                                allowed_states=[locations[s] for s in sorted(global_states.get(node, []))],
                                unknown_country=unknown, branch_length=node.branch_length))
        if node in parents:
            parent = parents[node]
            source = locations[assigned[parent]] if k else None
            branches.append(dict(parent_id=ids[parent], child_id=ids[node], source=source,
                                 target=state, changed=source != state, unknown_tip=unknown))
    buffer = StringIO()
    Phylo.write(tree, buffer, "newick")
    audit = dict(root_id=ids[tree.root], optimum_changes=optimum, tie_seed=seed,
                 tie_method="seeded local optimal-state choices; not uniform history sampling",
                 nodes=audit_nodes, branches=branches, newick=buffer.getvalue().strip())
    ranges = {pair: (int(minima[i]), int(maxima[i])) for i, pair in enumerate(pairs)
              if maxima[i] > 0}
    return counts, ranges, audit


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
    representative, ranges, audit = _reconstruct(tree, countries, locations)
    counts, ambiguity_counts, ambiguity = Counter(), Counter(), 0
    for root in roots:
        possible, ambiguous_pairs, root_ambiguity = _possible_pairs(tree, countries, locations, root)
        counts.update(possible)
        ambiguity_counts.update(ambiguous_pairs)
        ambiguity += root_ambiguity
    nodes = []
    for country in locations:
        subset = [r for r in rows if _country(r) == country]
        inputs = sum(r.get("origin") in _INPUT for r in subset)
        nodes.append(dict(country=country, count=len(subset), input_count=inputs,
                          context_count=len(subset) - inputs, is_input_country=bool(inputs)))
    return dict(
        method="weighted undirected changes in one coherent optimal parsimony history",
        reconstruction=audit, country_colors=country_palette(locations),
        optimum_changes=audit["optimum_changes"],
        interpretation="Representative counts describe one optimal history; exact pair ranges span "
        "all optimal histories on the fixed tree and need not be jointly attainable. "
        "Location changes are not demonstrated transmission. Root fractions "
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
        edges=[dict(source=a, target=b, representative_count=representative[(a, b)],
                    min_changes=ranges[(a, b)][0], max_changes=ranges[(a, b)][1],
                    roots_with_possible_change=count,
                    roots_tested=len(roots), root_fraction=count / len(roots),
                    uncertain=True, count_ambiguous=ranges[(a, b)][0] != ranges[(a, b)][1],
                    roots_with_ambiguous_change=ambiguity_counts[(a, b)])
               for (a, b) in sorted(ranges) for count in [counts[(a, b)]]],
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
        data["country_colors"] = dict(network["country_colors"])
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


def _weighted_layout(nodes, edges, seed=1729):
    """Deterministic NumPy spring layout; attraction uses representative counts."""
    import numpy as np

    names = [node["country"] for node in nodes]
    n = len(names)
    if not n:
        return {}
    if n == 1:
        return {names[0]: (0.0, 0.0)}
    index = {name: i for i, name in enumerate(names)}
    weights = np.zeros((n, n))
    for edge in edges:
        i, j = index[edge["source"]], index[edge["target"]]
        weights[i, j] = weights[j, i] = 0.65 * math.log1p(edge.get("representative_count", 0))
    positions = np.random.default_rng(seed).uniform(-1, 1, (n, 2))
    ideal = math.sqrt(1 / n)
    for iteration in range(160):
        delta = positions[:, None, :] - positions[None, :, :]
        distance = np.maximum(np.linalg.norm(delta, axis=2), 0.01)
        force = ideal**2 / distance**2 - weights * distance / ideal
        np.fill_diagonal(force, 0)
        displacement = np.sum(delta * force[:, :, None], axis=1)
        norm = np.maximum(np.linalg.norm(displacement, axis=1), 1e-9)
        temperature = 0.12 * (1 - iteration / 160)
        positions += displacement / norm[:, None] * np.minimum(norm, temperature)[:, None]
        positions -= positions.mean(axis=0)
    positions /= max(float(np.abs(positions).max()), 1e-9)
    # Give dense hubs room for node markers; preserve the force-layout geometry.
    for _ in range(60):
        delta = positions[:, None, :] - positions[None, :, :]
        distance = np.maximum(np.linalg.norm(delta, axis=2), 1e-9)
        overlap = np.maximum(0.18 - distance, 0)
        np.fill_diagonal(overlap, 0)
        positions += np.sum(delta / distance[:, :, None] * overlap[:, :, None], axis=1) * 0.25
    return {name: tuple(positions[i]) for i, name in enumerate(names)}


def _place_country_labels(fig, ax, nodes, positions):
    """Choose label offsets against measured text boxes and all node markers."""
    from matplotlib.transforms import Bbox
    from matplotlib.font_manager import FontProperties

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    scale = fig.dpi / 72
    font = FontProperties(size=9)
    obstacles = []
    for node in nodes:
        x, y = ax.transData.transform(positions[node["country"]])
        radius = math.sqrt(60 + 6 * node["count"]) / 2 * scale + 3
        obstacles.append(Bbox.from_extents(x - radius, y - radius, x + radius, y + radius))
    labels = []
    for node in sorted(nodes, key=lambda n: (-n["count"], n["country"])):
        name = node["country"]
        text = f"{name} ({node['count']})"
        width, height, _ = renderer.get_text_width_height_descent(text, font, ismath=False)
        height += 3
        px, py = ax.transData.transform(positions[name])
        radius = math.sqrt(60 + 6 * node["count"]) / 2 + 5
        candidates = []
        for extra in (0, 12, 26, 44, 66):
            for dx, dy in ((1, 1), (1, -1), (-1, 1), (-1, -1),
                           (1, 0), (-1, 0), (0, 1), (0, -1)):
                ox, oy = dx * (radius + extra), dy * (radius + extra)
                anchor_x, anchor_y = px + ox * scale, py + oy * scale
                left = anchor_x if dx >= 0 else anchor_x - width
                bottom = anchor_y - height / 2
                box = Bbox.from_bounds(left - 3, bottom - 2, width + 6, height + 4)
                intersections = [Bbox.intersection(box, other) for other in obstacles + labels]
                overlap = sum(b.width * b.height for b in intersections if b is not None)
                bounds = ax.get_window_extent()
                outside = max(bounds.x0 - box.x0, 0) + max(box.x1 - bounds.x1, 0)
                outside += max(bounds.y0 - box.y0, 0) + max(box.y1 - bounds.y1, 0)
                candidates.append((overlap + outside * 1000 + extra * 0.01,
                                   ox, oy, dx, box, extra))
        _, ox, oy, dx, box, extra = min(candidates, key=lambda candidate: candidate[0])
        labels.append(box)
        ax.annotate(text, positions[name], xytext=(ox, oy), textcoords="offset points",
                    ha="left" if dx >= 0 else "right", va="center", fontsize=9, zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.8),
                    arrowprops=dict(arrowstyle="-", color="#999999", linewidth=0.5)
                    if extra else None)


def draw_profile_network(network, path, focus_inputs=False, include_possible=False):
    """Plot weighted coherent history; exact numerical counts live in the audit/table."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    inputs = set(network["input_countries"])
    edges = [e for e in network["edges"] if (include_possible or e["representative_count"] > 0)
             and (not focus_inputs or e["source"] in inputs or e["target"] in inputs)]
    visible = inputs | {e[key] for e in edges for key in ("source", "target")}
    nodes = [n for n in network["nodes"] if not focus_inputs or n["country"] in visible]
    positions = _weighted_layout(nodes, edges)
    colours = network.get("country_colors") or country_palette(n["country"] for n in nodes)
    fig, ax = plt.subplots(figsize=(12, 10))
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.18, top=0.93)
    for edge in edges:
        a, b = positions[edge["source"]], positions[edge["target"]]
        count = edge["representative_count"]
        ambiguous = edge["min_changes"] != edge["max_changes"]
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#888888", alpha=0.6 if count else 0.25,
                linewidth=0.7 + 1.2 * math.sqrt(count),
                linestyle="--" if ambiguous else "-", zorder=1)
    for node in nodes:
        x, y = positions[node["country"]]
        ax.scatter(x, y, s=60 + 6 * node["count"], color=colours[node["country"]],
                   edgecolors="#222222" if node["is_input_country"] else "white",
                   linewidths=2 if node["is_input_country"] else 0.6, zorder=2)
    if nodes:
        xs, ys = zip(*positions.values())
        ax.set(xlim=(min(xs) - 0.23, max(xs) + 0.27),
               ylim=(min(ys) - 0.15, max(ys) + 0.15))
        _place_country_labels(fig, ax, nodes, positions)
    else:
        ax.text(0.5, 0.5, "No known country observations", transform=ax.transAxes, ha="center")
    handles = [Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#222222",
                      markeredgewidth=2, linestyle="", label="Contains input genomes"),
               Line2D([], [], color="#888888", label="Count fixed across optimal histories"),
               Line2D([], [], color="#888888", linestyle="--", label="Count varies across optimal histories")]
    for count in (1, 3, 10):
        handles.append(Line2D([], [], color="#888888", linewidth=0.7 + 1.2 * math.sqrt(count),
                              label=f"{count} reconstructed change" + ("s" if count != 1 else "")))
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.105),
               ncol=3, frameon=False, fontsize=9, columnspacing=2.0, handlelength=2.8)
    ax.set_title("Reconstructed country changes · " +
                 ("input-country focus" if focus_inputs else "all sampled countries"), fontsize=14, pad=18)
    fig.text(0.5, 0.035, "Node area: sampled genomes. Width: representative change count. "
             "Exact counts and ranges: accompanying table/CSV.\n"
             "One optimal history on a fixed tree; links are not demonstrated transmission. "
             f"Unknown tips: {network['unknown_country_count']} (wildcards).",
             ha="center", va="center", fontsize=9, linespacing=1.6)
    ax.axis("off")
    fig.savefig(path)
    plt.close(fig)
    return str(path)
