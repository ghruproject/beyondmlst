"""cgMLST LIN-prefix views over the shared biological-tree location network."""

import copy

from chronoclade.context_refinement import _lineage
from chronoclade.typing_scopes import compatible_typing as _compatible
from chronoclade.location_network.reconstruction import (
    build_location_network,
    is_input_sample,
    nearest_country_connections,
)


def build_profile_network(tree, rows, nearest_neighbours=None, subgroup_depth=5):
    """Build full-pool network and exact, version-compatible input LIN-prefix views."""
    network = build_location_network(tree, rows, nearest_neighbours, tree_basis="cgmlst")
    groups = {}
    for row in rows:
        lineage = _lineage(row, "cglin", subgroup_depth)
        if is_input_sample(row) and lineage is not None:
            groups.setdefault(lineage, []).append(row)
    views = []
    for index, (lineage, inputs) in enumerate(sorted(groups.items()), 1):
        subset = [
            row
            for row in rows
            if _lineage(row, "cglin", subgroup_depth) == lineage
            and _compatible(inputs[0], row, "cglin")
        ]
        ids = {r["sample_id"] for r in subset}
        pruned = copy.deepcopy(tree)
        for tip in list(pruned.get_terminals()):
            if tip.name not in ids:
                pruned.prune(tip)
        data = build_location_network(pruned, subset, nearest_neighbours, tree_basis="cgmlst")
        data["country_colors"] = dict(network["country_colors"])
        input_ids = set(data["input_ids"])
        subgroup_nearest = [
            item for item in nearest_neighbours or [] if item.get("query_id") in input_ids
        ]
        data["nearest_edges"] = nearest_country_connections(rows, subgroup_nearest)
        data["nearest_search_scope"] = "complete selected CG pool; subgroup input queries"
        prefix = ",".join(map(str, lineage[-1]))
        views.append(
            dict(
                id=f"lin_{index}",
                label=f"LIN level {subgroup_depth}: {prefix} ({len(inputs)} inputs)",
                scope=dict(
                    kind="cglin",
                    level=subgroup_depth,
                    prefix=prefix,
                    scheme=lineage[0],
                    version=lineage[1],
                ),
                **data,
            )
        )
    network["views"] = views
    return network
