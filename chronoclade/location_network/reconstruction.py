"""Directed country-state changes and exact optimal-history count ranges.

Ranges condition on a fixed tree. Root agreement is a separate sensitivity audit,
never transmission probability.
Nearest-country comparisons are observations and remain separate from reconstruction.
"""

from __future__ import annotations

import copy
import math
from collections import Counter
from typing import Literal

from Bio.Phylo.BaseTree import Tree

from .colours import country_palette

INPUT_ORIGINS = frozenset({"local", "query", "focal", "input"})


def is_input_sample(row):
    """Prepared sample roles take precedence over legacy source origins."""
    return (row.get("role") or row.get("origin")) in INPUT_ORIGINS


def _country(row):
    value = str(row.get("country") or "").strip()
    return None if value.casefold() in {"", "unknown", "none", "na", "n/a"} else value


def _roots(tree, rows, limit=5):
    """Seed with an input; prioritise new countries, then farthest tree tips."""
    by_id = {r["sample_id"]: r for r in rows}
    candidates = sorted(t.name for t in tree.get_terminals() if t.name in by_id)
    if not candidates:
        return []
    inputs = [name for name in candidates if is_input_sample(by_id[name])]
    positions = {name: index for index, name in enumerate(candidates)}
    selected = [inputs[0] if inputs else candidates[0]]
    while len(selected) < min(limit, len(candidates)):
        seen = {_country(by_id[name]) for name in selected} - {None}
        remaining = [name for name in candidates if name not in selected]
        novel = [
            name
            for name in remaining
            if _country(by_id[name]) not in seen and _country(by_id[name]) is not None
        ]
        choices = novel or remaining
        # Stable tie-breaking by accession; actual selection is tree-distance based.
        selected.append(
            max(
                choices,
                key=lambda name: (
                    min(tree.distance(name, root) for root in selected),
                    -positions[name],
                ),
            )
        )
    return selected


def _possible_pairs(tree, countries, locations, root):
    rooted = copy.deepcopy(tree)
    if len(rooted.get_terminals()) > 1:
        rooted.root_with_outgroup(root)
    costs = {}
    for node in rooted.find_clades(order="postorder"):
        if node.is_terminal():
            observed = countries.get(node.name)
            costs[node] = {
                state: 0 if observed in {state, None} else math.inf for state in locations
            }
        else:
            costs[node] = {
                state: sum(
                    min(costs[child][target] + (state != target) for target in locations)
                    for child in node.clades
                )
                for state in locations
            }
    optimum = min(costs[rooted.root].values())
    states = {rooted.root: {s for s in locations if costs[rooted.root][s] == optimum}}
    possible, ambiguous_pairs, ambiguous = set(), set(), 0
    for parent in rooted.find_clades(order="preorder"):
        for child in parent.clades:
            allowed = set()
            for source in states[parent]:
                best = min(costs[child][target] + (source != target) for target in locations)
                allowed.update(
                    (source, target)
                    for target in locations
                    if costs[child][target] + (source != target) == best
                )
            states[child] = {target for _, target in allowed}
            ambiguous += len(allowed) > 1
            pairs = {tuple(sorted((a, b))) for a, b in allowed if a != b}
            possible.update(pairs)
            if len(allowed) > 1:
                ambiguous_pairs.update(pairs)
    return possible, ambiguous_pairs, ambiguous


def _reconstruct(tree, countries, locations, seed=1729, directed=False):
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
    pairs = (
        [(a, b) for a in locations for b in locations if a != b]
        if directed
        else [(a, b) for i, a in enumerate(locations) for b in locations[i + 1 :]]
    )
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
                        pair = (locations[source], locations[target])
                        col = pair_index[pair if directed else tuple(sorted(pair))]
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
                global_states[child] = {
                    int(t) for s in global_states[node] for t in allowed[child][s]
                }
                if source != target:
                    pair = (locations[source], locations[target])
                    counts[pair if directed else tuple(sorted(pair))] += 1
    audit_nodes, branches = [], []
    for node in nodes:
        state = locations[assigned[node]] if k else None
        observed = countries.get(node.name) if node.is_terminal() else None
        unknown = node.is_terminal() and observed is None
        audit_nodes.append(
            dict(
                id=ids[node],
                parent_id=ids[parents[node]] if node in parents else None,
                children_ids=[ids[c] for c in node.clades],
                name=node.name,
                is_tip=node.is_terminal(),
                state=None if unknown else state,
                inferred_state=state,
                observed_country=observed,
                allowed_states=[locations[s] for s in sorted(global_states.get(node, []))],
                unknown_country=unknown,
                branch_length=node.branch_length,
            )
        )
        if node in parents:
            parent = parents[node]
            source = locations[assigned[parent]] if k else None
            branches.append(
                dict(
                    parent_id=ids[parent],
                    child_id=ids[node],
                    source=source,
                    target=state,
                    changed=source != state,
                    unknown_tip=unknown,
                )
            )
    buffer = StringIO()
    Phylo.write(tree, buffer, "newick")
    audit = dict(
        root_id=ids[tree.root],
        optimum_changes=optimum,
        tie_seed=seed,
        tie_method="seeded local optimal-state choices; not uniform history sampling",
        nodes=audit_nodes,
        branches=branches,
        newick=buffer.getvalue().strip(),
    )
    ranges = {
        pair: (int(minima[i]), int(maxima[i])) for i, pair in enumerate(pairs) if maxima[i] > 0
    }
    return counts, ranges, audit


def nearest_country_connections(rows, comparisons):
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
        row = dict(
            source=source,
            target=target,
            query_ids=qs,
            context_ids=cs,
            query_count=len(qs),
            context_count=len(cs),
            comparison_count=len(items),
        )
        for field, alias in (
            ("allele_differences", "allele_mismatches"),
            ("shared_called_loci", "shared_loci"),
        ):
            values = [i[field] for i in items if i.get(field) is not None]
            row[alias + "_min"] = min(values) if values else None
            row[alias + "_max"] = max(values) if values else None
        result.append(row)
    return result


def _network_metrics(locations, edges):
    """Unweighted directed graph metrics plus weighted representative counts."""
    representative = [edge for edge in edges if edge["representative_count"] > 0]
    vertices = list(locations)
    outgoing = {country: set() for country in vertices}
    incoming = {country: set() for country in vertices}
    weights_out = {country: 0 for country in vertices}
    weights_in = {country: 0 for country in vertices}
    for edge in representative:
        source, target = edge["source"], edge["target"]
        outgoing[source].add(target)
        incoming[target].add(source)
        weights_out[source] += edge["representative_count"]
        weights_in[target] += edge["representative_count"]

    # Brandes' algorithm for directed, unweighted betweenness. Like igraph's
    # default, values are unnormalised and paths follow edge direction.
    betweenness = {country: 0.0 for country in vertices}
    for source in vertices:
        stack, predecessors = [], {country: [] for country in vertices}
        paths = {country: 0 for country in vertices}
        distance = {country: -1 for country in vertices}
        paths[source], distance[source] = 1, 0
        queue = [source]
        for vertex in queue:
            stack.append(vertex)
            for neighbour in sorted(outgoing[vertex]):
                if distance[neighbour] < 0:
                    queue.append(neighbour)
                    distance[neighbour] = distance[vertex] + 1
                if distance[neighbour] == distance[vertex] + 1:
                    paths[neighbour] += paths[vertex]
                    predecessors[neighbour].append(vertex)
        dependency = {country: 0.0 for country in vertices}
        while stack:
            vertex = stack.pop()
            for parent in predecessors[vertex]:
                dependency[parent] += (paths[parent] / paths[vertex]) * (1 + dependency[vertex])
            if vertex != source:
                betweenness[vertex] += dependency[vertex]

    # igraph closeness(mode="all"): traverse the undirected adjacency and use
    # the reciprocal of the sum of reachable distances, without normalization.
    adjacency = {country: outgoing[country] | incoming[country] for country in vertices}
    closeness = {}
    for source in vertices:
        distances = {source: 0}
        queue = [source]
        for vertex in queue:
            for neighbour in sorted(adjacency[vertex]):
                if neighbour not in distances:
                    distances[neighbour] = distances[vertex] + 1
                    queue.append(neighbour)
        distance_sum = sum(distances.values())
        closeness[source] = 1 / distance_sum if distance_sum else None

    result = []
    for country in locations:
        in_degree, out_degree = len(incoming[country]), len(outgoing[country])
        total_degree = in_degree + out_degree
        in_weight, out_weight = weights_in[country], weights_out[country]
        result.append(
            dict(
                country=country,
                in_degree=in_degree,
                out_degree=out_degree,
                degree=total_degree,
                betweenness=betweenness[country],
                closeness=closeness[country],
                in_changes=in_weight,
                out_changes=out_weight,
                source_hub_ratio=out_degree / total_degree if total_degree else None,
            )
        )
    return result


def build_location_network(
    tree: Tree,
    rows,
    nearest_neighbours=None,
    *,
    tree_basis: Literal["cgmlst", "sequence", "provided_phylogeny"],
):
    """Reconstruct locations on an explicitly supplied biological tree.

    A protein embedding or an embedding-distance NJ diagram is not a biological
    phylogeny. ESM2 workflows can use this component only with an independently
    supplied biological tree, declared as ``provided_phylogeny`` or ``sequence``.
    Tree provenance is the caller's responsibility; it cannot be inferred from
    a Newick topology alone. No tree is constructed by this component.
    """
    if tree_basis not in {"cgmlst", "sequence", "provided_phylogeny"}:
        raise ValueError(
            "Location reconstruction requires a biological tree, not embedding distances"
        )
    if not isinstance(tree, Tree):
        raise TypeError("Location reconstruction requires a Bio.Phylo biological tree")
    tip_ids = {tip.name for tip in tree.get_terminals()}
    rows = [r for r in rows if r["sample_id"] in tip_ids]
    countries = {r["sample_id"]: _country(r) for r in rows}
    locations = sorted(set(countries.values()) - {None})
    roots = _roots(tree, rows) if locations else []
    representative, ranges, audit = _reconstruct(tree, countries, locations)
    directed_counts, directed_ranges, directed_audit = _reconstruct(
        tree, countries, locations, directed=True
    )
    assert audit == directed_audit  # Both objectives use precisely the same seeded history.
    assert sum(directed_counts.values()) == audit["optimum_changes"]
    directed_edges = [
        dict(
            source=a,
            target=b,
            representative_count=directed_counts[(a, b)],
            min_changes=low,
            max_changes=high,
            uncertain=True,
            count_ambiguous=low != high,
        )
        for (a, b), (low, high) in sorted(directed_ranges.items())
    ]
    counts, ambiguity_counts, ambiguity = Counter(), Counter(), 0
    for root in roots:
        possible, ambiguous_pairs, root_ambiguity = _possible_pairs(
            tree, countries, locations, root
        )
        counts.update(possible)
        ambiguity_counts.update(ambiguous_pairs)
        ambiguity += root_ambiguity
    nodes = []
    for country in locations:
        subset = [r for r in rows if _country(r) == country]
        inputs = sum(is_input_sample(r) for r in subset)
        nodes.append(
            dict(
                country=country,
                count=len(subset),
                input_count=inputs,
                context_count=len(subset) - inputs,
                is_input_country=bool(inputs),
            )
        )
    return dict(
        tree_basis=tree_basis,
        method="weighted directed parent-to-child country-state changes in one optimal parsimony history",
        directed_edges=directed_edges,
        network_metrics=_network_metrics(locations, directed_edges),
        direction_scope=(
            "original rooted NJ tree" if tree_basis == "cgmlst" else "original rooted tree"
        )
        + "; directions can change under rerooting",
        reconstruction=audit,
        country_colors=country_palette(locations),
        optimum_changes=audit["optimum_changes"],
        interpretation="Representative counts describe one optimal history; exact pair ranges span "
        "all optimal histories on the fixed rooted tree and need not be jointly attainable. "
        "Arrows follow parent-to-child country-state changes; directions depend on the root. "
        "Location changes are not demonstrated transmission. Root fractions "
        "measure sensitivity on a fixed tree, not probability or confidence; sampling and tree "
        "uncertainty are not estimated.",
        root_selection="input seed, then country diversity and farthest tree tips",
        uncertainty_scope="fixed tree; no topology or sampling uncertainty assessment",
        tested_roots=roots,
        root_count=len(roots),
        ambiguous_edge_reconstructions=ambiguity,
        nodes=nodes,
        unknown_country_count=sum(value is None for value in countries.values()),
        input_countries=[n["country"] for n in nodes if n["is_input_country"]],
        input_ids=sorted(r["sample_id"] for r in rows if is_input_sample(r)),
        sample_count=len(rows),
        edges=[
            dict(
                source=a,
                target=b,
                representative_count=representative[(a, b)],
                min_changes=ranges[(a, b)][0],
                max_changes=ranges[(a, b)][1],
                roots_with_possible_change=count,
                roots_tested=len(roots),
                root_fraction=count / len(roots),
                uncertain=True,
                count_ambiguous=ranges[(a, b)][0] != ranges[(a, b)][1],
                roots_with_ambiguous_change=ambiguity_counts[(a, b)],
            )
            for (a, b) in sorted(ranges)
            for count in [counts[(a, b)]]
        ],
        nearest_edges=nearest_country_connections(rows, nearest_neighbours),
    )
