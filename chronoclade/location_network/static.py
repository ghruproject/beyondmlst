"""Weighted, reproducible network figures from one coherent history."""

import math

from .colours import country_palette


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
        weights[i, j] += edge.get("representative_count", 0)
        weights[j, i] += edge.get("representative_count", 0)
    weights = 0.65 * np.log1p(weights)
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
    font_size = 11 if len(nodes) <= 16 else 9
    font = FontProperties(size=font_size)
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
            for dx, dy in ((1, 1), (1, -1), (-1, 1), (-1, -1), (1, 0), (-1, 0), (0, 1), (0, -1)):
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
                candidates.append((overlap + outside * 1000 + extra * 0.01, ox, oy, dx, box, extra))
        _, ox, oy, dx, box, extra = min(candidates, key=lambda candidate: candidate[0])
        labels.append(box)
        ax.annotate(
            text,
            positions[name],
            xytext=(ox, oy),
            textcoords="offset points",
            ha="left" if dx >= 0 else "right",
            va="center",
            fontsize=font_size,
            zorder=4,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.8),
            arrowprops=dict(arrowstyle="-", color="#999999", linewidth=0.5) if extra else None,
        )


def draw_location_network(network, path, focus_inputs=False, include_possible=False):
    """Plot weighted coherent history; exact numerical counts live in the audit/table."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch

    directed = "directed_edges" in network
    inputs = set(network["input_countries"])
    edges = [
        e
        for e in network["directed_edges" if directed else "edges"]
        if (include_possible or e["representative_count"] > 0)
        and (not focus_inputs or e["source"] in inputs or e["target"] in inputs)
    ]
    visible = inputs | {e[key] for e in edges for key in ("source", "target")}
    nodes = [n for n in network["nodes"] if not focus_inputs or n["country"] in visible]
    positions = _weighted_layout(nodes, edges)
    colours = network.get("country_colors") or country_palette(n["country"] for n in nodes)
    fig, ax = plt.subplots(figsize=(12, 10))
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.18, top=0.93)
    pairs = {(e["source"], e["target"]) for e in edges}
    node_counts = {n["country"]: n["count"] for n in nodes}
    for edge in edges:
        a, b = positions[edge["source"]], positions[edge["target"]]
        count = edge["representative_count"]
        ambiguous = edge["min_changes"] != edge["max_changes"]
        style = dict(
            color="#888888",
            alpha=0.65 if count else 0.25,
            linewidth=0.7 + 1.2 * math.sqrt(count),
            linestyle="--" if ambiguous else "-",
            zorder=1,
        )
        if directed:
            reciprocal = (edge["target"], edge["source"]) in pairs
            ax.add_patch(
                FancyArrowPatch(
                    a,
                    b,
                    arrowstyle="-|>",
                    mutation_scale=13 + math.sqrt(count),
                    connectionstyle=f"arc3,rad={0.16 if reciprocal else 0.035}",
                    shrinkA=math.sqrt(60 + 6 * node_counts[edge["source"]]) / 2 + 2,
                    shrinkB=math.sqrt(60 + 6 * node_counts[edge["target"]]) / 2 + 2,
                    **style,
                )
            )
        else:
            ax.plot([a[0], b[0]], [a[1], b[1]], **style)
    for node in nodes:
        x, y = positions[node["country"]]
        ax.scatter(
            x,
            y,
            s=60 + 6 * node["count"],
            color=colours[node["country"]],
            edgecolors="#222222" if node["is_input_country"] else "white",
            linewidths=2 if node["is_input_country"] else 0.6,
            zorder=2,
        )
    if nodes:
        xs, ys = zip(*positions.values())
        ax.set(xlim=(min(xs) - 0.23, max(xs) + 0.27), ylim=(min(ys) - 0.15, max(ys) + 0.15))
        _place_country_labels(fig, ax, nodes, positions)
    else:
        ax.text(0.5, 0.5, "No known country observations", transform=ax.transAxes, ha="center")
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            markerfacecolor="white",
            markeredgecolor="#222222",
            markeredgewidth=2,
            linestyle="",
            label="Contains input genomes",
        ),
        Line2D([], [], color="#888888", label="Count fixed across optimal histories"),
        Line2D(
            [], [], color="#888888", linestyle="--", label="Count varies across optimal histories"
        ),
    ]
    for count in (1, 3, 10):
        handles.append(
            Line2D(
                [],
                [],
                color="#888888",
                linewidth=0.7 + 1.2 * math.sqrt(count),
                label=f"{count} reconstructed change" + ("s" if count != 1 else ""),
            )
        )
    fig.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.105),
        ncol=3,
        frameon=False,
        fontsize=9,
        columnspacing=2.0,
        handlelength=2.8,
    )
    ax.set_title(
        ("Directed country-state changes · " if directed else "Reconstructed country changes · ")
        + ("input-country focus" if focus_inputs else "all sampled countries"),
        fontsize=14,
        pad=18,
    )
    direction_note = (
        "Arrows: parent → child states on the rooted tree; root-dependent, not proven transmission. "
        if directed
        else "One optimal history on a fixed tree; links are not demonstrated transmission. "
    )
    fig.text(
        0.5,
        0.035,
        "Node area: sampled genomes. Width: representative change count. "
        "Exact counts and ranges: accompanying table/CSV.\n"
        + direction_note
        + f"Unknown tips: {network['unknown_country_count']} (wildcards).",
        ha="center",
        va="center",
        fontsize=9,
        linespacing=1.6,
    )
    ax.axis("off")
    fig.savefig(path)
    plt.close(fig)
    return str(path)
