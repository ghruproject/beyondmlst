"""Country annotations from the same coherent history used by the fast network."""

from pathlib import Path

from chronoclade.sample_labels import sample_labels


def draw_country_tree(network, records, path, *, display_ids=None, nearest_ids=()):
    """Retain original ancestors when hiding tips; never reconstruct a display subset."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from chronoclade.profile_analysis import _date_interval

    reconstruction = network.get("reconstruction", {})
    nodes = {node["id"]: node for node in reconstruction.get("nodes", [])}
    path = Path(path)
    if not nodes:
        return None
    by_id = {row["sample_id"]: row for row in records}
    labels = sample_labels(records)
    tips = [node for node in nodes.values() if node["is_tip"]]
    shown = set(display_ids) if display_ids else {node["name"] for node in tips}
    shown &= {node["name"] for node in tips}
    nearest_ids = set(nearest_ids)
    root_id = reconstruction["root_id"]
    included = set()
    for node in tips:
        if node["name"] in shown:
            current = node["id"]
            while current is not None and current not in included:
                included.add(current)
                current = nodes[current]["parent_id"]
    if not included:
        return None
    xs, ys = {}, {}
    terminal_order = []

    def position(ident, depth):
        node = nodes[ident]
        xs[ident] = depth
        children = [child for child in node["children_ids"] if child in included]
        if not children:
            ys[ident] = len(terminal_order)
            terminal_order.append(ident)
        else:
            for child in children:
                position(child, depth + (nodes[child].get("branch_length") or 0))
            ys[ident] = sum(ys[child] for child in children) / len(children)

    position(root_id, 0)
    palette = network.get("country_colors", {})
    fig, ax = plt.subplots(figsize=(14, max(5, len(shown) * 0.25)))
    for ident in included:
        node = nodes[ident]
        children = [child for child in node["children_ids"] if child in included]
        colour = palette.get(node.get("state"), "#999999")
        if children:
            ax.plot([xs[ident]] * 2, [min(ys[c] for c in children), max(ys[c] for c in children)],
                    color=colour, linewidth=0.8)
        parent = node["parent_id"]
        if parent in included:
            ax.plot([xs[parent], xs[ident]], [ys[ident]] * 2, color=colour, linewidth=1,
                    linestyle="--" if (node.get("branch_length") or 0) < 0 else "-")
        ambiguous = len(node.get("allowed_states", [])) > 1
        if not node["is_tip"]:
            ax.plot(xs[ident], ys[ident], "o", markersize=3,
                    markerfacecolor="white" if ambiguous else colour, markeredgecolor=colour)
    span = max(xs.values()) - min(xs.values()) or 1
    label_x = max(xs.values()) + span * 0.025
    for ident in terminal_order:
        node = nodes[ident]
        row = by_id.get(node["name"], {})
        country = str(row.get("country") or "Unknown")
        colour = palette.get(country, "#999999")
        is_input = row.get("origin") in {"local", "query", "focal", "input"}
        ax.plot(xs[ident], ys[ident], "o", markersize=5 if is_input else 3,
                markerfacecolor=colour, markeredgecolor="#111111" if is_input else colour,
                markeredgewidth=1.3 if is_input else 0.5, zorder=3)
        if node["name"] in nearest_ids:
            ax.plot(xs[ident], ys[ident], marker="*", markersize=9, markerfacecolor="none",
                    markeredgecolor="#e66101", linestyle="None", zorder=4)
        date = str(row.get("collection_date") or row.get("collection_year") or "Date unknown")
        if date != "Date unknown" and _date_interval(row) is None:
            date += " (excluded from date analysis)"
        text = f"{labels.get(node['name'], node['name'])} | {country} | {date}"
        ax.plot([xs[ident], label_x], [ys[ident]] * 2, color="#dddddd", linewidth=0.4)
        ax.text(label_x, ys[ident], text, va="center", fontsize=8,
                fontweight="bold" if is_input else "normal", color=colour)
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor="#bbbbbb",
                     markeredgecolor="#111111", label="Outlined tips: input genomes"),
               Line2D([], [], marker="o", color="#777777", markerfacecolor="white",
                      linestyle="None", label="Hollow ancestors: multiple optimal countries"),
               Line2D([], [], marker="*", color="#e66101", markerfacecolor="none",
                      linestyle="None", label="Nearest public relatives (ties retained)")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False,
              fontsize=9, ncol=1)
    ax.set_title("Country-coloured NJ tree · same reconstruction as the weighted network", pad=70)
    ax.set_xlabel("Fraction of mismatching callable cgMLST loci; dashed branches have negative length")
    ax.set_yticks([])
    ax.set_ylim(len(terminal_order), -1)
    ax.set_xlim(min(xs.values()) - span * 0.03, max(xs.values()) + span * 1.15)
    ax.text(0, -0.04, f"Showing {len(shown)} of {len(tips)} profiles. Ancestral colours show one optimal "
            "history; hollow nodes retain ambiguity. This tree is not dated.", transform=ax.transAxes,
            fontsize=9)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return str(path)
