"""Auditable country changes in TreeTime's marginal ancestral reconstruction."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import textwrap
from collections import Counter, defaultdict
from html import escape
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError

METHOD = (
    "TreeTime marginal maximum-likelihood ancestral location reconstruction. Arrows count "
    "country-changing parent-to-child branches among the analysed genomes. They are not direct "
    "transmission links, dated journeys, or StrainHub parsimony results."
)
CAVEAT = (
    "Direction depends on the chosen root and sampled genomes. Tied ancestral assignments can "
    "change or reverse arrows. Endpoint marginal confidence is not joint support for an edge; "
    "this diagram does not quantify uncertainty over alternative trees, roots or sampling. "
    "Labels use locations as supplied: this is a country network only when those labels are countries."
)
UNKNOWN = {
    "",
    "?",
    "unknown",
    "unknown country",
    "public_context",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "missing",
}
BRANCH_FIELDS = [
    "parent",
    "child",
    "parent_country",
    "child_country",
    "parent_marginal_confidence",
    "child_marginal_confidence",
    "parent_candidates",
    "child_candidates",
    "parent_tied",
    "child_tied",
    "changes_country",
    "both_endpoints_high_confidence",
    "included",
    "reason",
]
EDGE_FIELDS = ["source", "target", "branches", "high_endpoint_branches", "uncertain_branches"]
NODE_FIELDS = ["country", "observed_tips"]


def _known(value: str) -> bool:
    return value.strip().casefold() not in UNKNOWN


def _csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _mapping(path: Path) -> dict[str, str]:
    mapping = {}
    active = False
    for line in path.read_text().splitlines():
        if line.strip() == "Character to attribute mapping:":
            active = True
            continue
        if active and not line.strip():
            break
        if active:
            match = re.fullmatch(r"\s*(\S+):\s*(.+?)\s*", line)
            if not match or match[1] in mapping:
                raise ValueError("Malformed or duplicate location mapping")
            mapping[match[1]] = match[2]
    if not mapping:
        raise ValueError("Location character mapping is unavailable")
    return mapping


def _confidences(path: Path, mapping: dict) -> dict:
    values = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, skipinitialspace=True)
        header = next(reader, [])
        if not header:
            raise ValueError("Location confidence table is empty")
        codes = [x.strip() for x in header[1:]]
        if header[0].strip() != "#name" or not codes or len(codes) != len(set(codes)):
            raise ValueError("Malformed location confidence header")
        if any(code not in mapping for code in codes):
            raise ValueError("Unmapped confidence character")
        for row in reader:
            if len(row) != len(header) or not row[0].strip() or row[0].strip() in values:
                raise ValueError("Malformed or duplicate confidence row")
            probabilities = [float(x) for x in row[1:]]
            if any(not math.isfinite(x) or x < 0 or x > 1 for x in probabilities):
                raise ValueError("Invalid location probability")
            if not math.isclose(sum(probabilities), 1, abs_tol=1e-6):
                raise ValueError("Location probabilities do not sum to one")
            states = {mapping[code]: probability for code, probability in zip(codes, probabilities)}
            if len(states) != len(codes):
                raise ValueError("Duplicate country labels in confidence mapping")
            maximum = max(probabilities)
            candidates = sorted(
                state
                for state, value in states.items()
                if math.isclose(value, maximum, rel_tol=0, abs_tol=1e-10)
            )
            values[row[0].strip()] = {"probabilities": states, "candidates": candidates}
    return values


def _state(clade, confidences: dict) -> dict:
    name = clade.name
    match = re.search(r'(?:^|[,\[&])location="([^"]*)"', clade.comment or "")
    if not name or not match or name not in confidences:
        raise ValueError("Tree node is missing its name, location annotation or confidence row")
    country = match[1]
    confidence = confidences[name]
    if country not in confidence["probabilities"]:
        # TreeTime's missing-state annotation is excluded, not assigned to a real country.
        if not _known(country):
            return {"country": country, "confidence": None, "candidates": [], "tied": False}
        raise ValueError(f"Unmapped annotated country for node {name}")
    if country not in confidence["candidates"]:
        raise ValueError(f"Annotated country disagrees with marginal maximum at node {name}")
    return {
        "country": country,
        "confidence": confidence["probabilities"][country],
        "candidates": confidence["candidates"],
        "tied": len(confidence["candidates"]) > 1,
    }


def _figure(directory: Path, result: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch

    nodes = result["nodes"]
    fig, ax = plt.subplots(figsize=(9, max(6.5, len(nodes) * 0.32)))
    ax.set_axis_off()
    ax.set_xlim(-1.85, 1.85)
    ax.set_ylim(-1.5, 1.5)
    ax.set_title("Location changes on the analysed tree", fontsize=16, pad=16)
    if result["status"] == "unavailable":
        ax.text(
            0,
            0,
            "Country reconstruction unavailable\n" + result["reason"],
            ha="center",
            va="center",
            wrap=True,
            fontsize=12,
        )
    elif not nodes:
        ax.text(0, 0, "No known country labels", ha="center", fontsize=13)
    else:
        positions = {
            node["country"]: (
                math.cos(2 * math.pi * i / len(nodes)),
                math.sin(2 * math.pi * i / len(nodes)),
            )
            for i, node in enumerate(nodes)
        }
        label_positions = {}
        if len(nodes) > 12:
            for side in (-1, 1):
                group = sorted(
                    (
                        country
                        for country, (x, y) in positions.items()
                        if (1 if x >= 0 else -1) == side
                    ),
                    key=lambda country: positions[country][1],
                )
                for index, country in enumerate(group):
                    label_positions[country] = (
                        side * 1.25,
                        -1.3 + 2.6 * index / max(1, len(group) - 1),
                    )
        for edge in result["edges"]:
            start, end = positions[edge["source"]], positions[edge["target"]]
            uncertain = edge["uncertain_branches"] > 0
            patch = FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=16,
                connectionstyle="arc3,rad=0.18",
                shrinkA=25,
                shrinkB=25,
                linewidth=1.5 + math.log1p(edge["branches"]),
                linestyle="--" if uncertain else "-",
                color="#7a629b" if uncertain else "#3f285d",
            )
            ax.add_patch(patch)
            midpoint = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
            ax.text(
                midpoint[0],
                midpoint[1],
                f"{edge['branches']} {'branch' if edge['branches'] == 1 else 'branches'}\n"
                f"{edge['uncertain_branches']} uncertain",
                fontsize=9,
                ha="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
            )
        for node in nodes:
            x, y = positions[node["country"]]
            ax.scatter(
                [x],
                [y],
                s=450 + 70 * node["observed_tips"],
                color="#e5ddef",
                edgecolor="#594074",
                zorder=3,
            )
            label_x, label_y = label_positions.get(node["country"], (x * 1.12, y * 1.22))
            if label_positions:
                ax.plot([x, label_x], [y, label_y], color="#c4b9d0", linewidth=0.6, zorder=1)
            label_country = textwrap.fill(node["country"], width=24)
            ax.text(
                label_x,
                label_y,
                f"{label_country}\n{node['observed_tips']} observed "
                f"{'genome' if node['observed_tips'] == 1 else 'genomes'}",
                ha="left" if label_x > 0.3 else "right" if label_x < -0.3 else "center",
                va="center",
                fontsize=10,
            )
        if not result["edges"]:
            ax.text(0, 0, "No reconstructed country changes", ha="center", fontsize=11)
    fig.text(
        0.5,
        0.025,
        "Solid: confident location assignments at both ends.\n"
        "Dashed: one or both endpoints uncertain. Arrows are not transmission.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(directory / "country_network.svg")
    fig.savefig(directory / "country_network.png", dpi=160)
    plt.close(fig)


def build_marginal_network(
    directory: Path, *, temporal_supported: bool = False, confidence_threshold: float = 0.9
) -> dict:
    """Write a network and complete branch audit; missing/invalid inputs give unavailable output."""
    directory = Path(directory)
    if not math.isfinite(confidence_threshold) or not 0 <= confidence_threshold <= 1:
        raise ValueError("Confidence threshold must be between zero and one")
    result = {
        "status": "unavailable",
        "reason": "",
        "method": METHOD,
        "caveat": CAVEAT,
        "temporal_supported": bool(temporal_supported),
        "confidence_threshold": confidence_threshold,
        "nodes": [],
        "edges": [],
        "branches": [],
        "excluded_tips": [],
        "unknown_tips": 0,
        "node_assignments": [],
        "input_sha256": {},
    }
    paths = {
        "metadata": directory / "metadata.csv",
        "tree": directory / "location/annotated_tree.nexus",
        "confidence": directory / "location/confidence.csv",
        "mapping": directory / "location/GTR.txt",
    }
    try:
        missing = [key for key, path in paths.items() if not path.is_file()]
        if "metadata" in missing:
            raise ValueError("Missing input: " + ", ".join(missing))
        result["input_sha256"] = {
            key: hashlib.sha256(path.read_bytes()).hexdigest()
            for key, path in paths.items()
            if path.is_file()
        }
        with paths["metadata"].open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not {"sample_id", "location"}.issubset(reader.fieldnames or []):
                raise ValueError("Metadata requires sample_id and location")
            metadata = {}
            for row in reader:
                if not isinstance(row["sample_id"], str) or not isinstance(row["location"], str):
                    raise ValueError("Malformed metadata row: sample_id and location are required")
                name = row["sample_id"].strip()
                if not name or name in metadata:
                    raise ValueError("Missing or duplicate metadata sample ID")
                metadata[name] = row["location"].strip()
        if len({value for value in metadata.values() if _known(value)}) < 2:
            raise ValueError(
                "Fewer than two known locations; no between-location reconstruction was performed"
            )
        if missing:
            raise ValueError("Missing input: " + ", ".join(missing))
        tree_text = paths["tree"].read_text()
        match = re.search(r"(?ims)^\s*tree\s+[^=]+=(.+?;)\s*$", tree_text)
        if not match:
            raise ValueError("Missing or malformed annotated tree")
        tree = Phylo.read(io.StringIO(match[1]), "newick")
        names = [node.name for node in tree.find_clades()]
        if any(not name or not name.strip() for name in names):
            raise ValueError("Missing tree node ID")
        if len(names) != len(set(names)):
            raise ValueError("Duplicate tree node IDs")
        if any(
            node.branch_length is not None
            and (not math.isfinite(node.branch_length) or node.branch_length < 0)
            for node in tree.find_clades()
        ):
            raise ValueError("Tree branch lengths must be finite and nonnegative")
        confidences = _confidences(paths["confidence"], _mapping(paths["mapping"]))
        states = {node.name: _state(node, confidences) for node in tree.find_clades()}
        result["node_assignments"] = [
            {
                "node_id": node.name,
                "country": country,
                "marginal_probability": probability,
                "maximum_candidate": country in confidences[node.name]["candidates"],
                "annotated_representative": country == states[node.name]["country"],
                "is_tip": node.is_terminal(),
            }
            for node in tree.find_clades()
            for country, probability in sorted(confidences[node.name]["probabilities"].items())
        ]
        excluded = {}
        observed = Counter()
        for tip in tree.get_terminals():
            if tip.name not in metadata:
                excluded[tip.name] = "No exact tip metadata match"
                result["excluded_tips"].append(
                    {"sample_id": tip.name, "reason": "No exact metadata match"}
                )
            elif not _known(metadata[tip.name]):
                result["unknown_tips"] += 1
                excluded[tip.name] = "Missing observed location"
                result["excluded_tips"].append(
                    {"sample_id": tip.name, "reason": "Missing observed location"}
                )
            elif metadata[tip.name] != states[tip.name]["country"]:
                raise ValueError(f"Tip annotation disagrees with metadata for {tip.name}")
            else:
                observed[metadata[tip.name]] += 1
        for name in sorted(set(metadata) - {tip.name for tip in tree.get_terminals()}):
            result["excluded_tips"].append(
                {"sample_id": name, "reason": "Absent from reconstruction tree"}
            )
        counts = defaultdict(Counter)
        for parent in tree.find_clades():
            for child in parent.clades:
                a, b = states[parent.name], states[child.name]
                reason = ""
                if child.name in excluded:
                    reason = excluded[child.name]
                elif not _known(a["country"]) or not _known(b["country"]):
                    reason = "Unknown endpoint country"
                high = (
                    not a["tied"]
                    and not b["tied"]
                    and a["confidence"] is not None
                    and b["confidence"] is not None
                    and a["confidence"] >= confidence_threshold
                    and b["confidence"] >= confidence_threshold
                )
                changes = a["country"] != b["country"]
                branch = {
                    "parent": parent.name,
                    "child": child.name,
                    "parent_country": a["country"],
                    "child_country": b["country"],
                    "parent_marginal_confidence": a["confidence"],
                    "child_marginal_confidence": b["confidence"],
                    "parent_candidates": " | ".join(a["candidates"]),
                    "child_candidates": " | ".join(b["candidates"]),
                    "parent_tied": a["tied"],
                    "child_tied": b["tied"],
                    "changes_country": changes,
                    "both_endpoints_high_confidence": high,
                    "included": not reason,
                    "reason": reason,
                }
                result["branches"].append(branch)
                if changes and not reason:
                    count = counts[(a["country"], b["country"])]
                    count["branches"] += 1
                    count["high_endpoint_branches" if high else "uncertain_branches"] += 1
        result["edges"] = [
            {
                "source": a,
                "target": b,
                "branches": c["branches"],
                "high_endpoint_branches": c["high_endpoint_branches"],
                "uncertain_branches": c["uncertain_branches"],
            }
            for (a, b), c in sorted(counts.items())
        ]
        countries = set(observed) | {state for pair in counts for state in pair}
        result["nodes"] = [
            {"country": country, "observed_tips": observed[country]}
            for country in sorted(countries)
        ]
        result["status"] = "partial" if result["excluded_tips"] else "available"
        result["changing_branches"] = sum(edge["branches"] for edge in result["edges"])
        result["uncertain_changing_branches"] = sum(
            edge["uncertain_branches"] for edge in result["edges"]
        )
    except (ValueError, OSError, StopIteration, KeyError, NewickError) as error:
        result.update(
            status="unavailable",
            reason=str(error),
            nodes=[],
            edges=[],
            branches=[],
            node_assignments=[],
        )
    result["unavailable_reason"] = result["reason"] if result["status"] == "unavailable" else ""
    directory.mkdir(parents=True, exist_ok=True)
    _csv(directory / "country_network_nodes.csv", NODE_FIELDS, result["nodes"])
    _csv(directory / "country_network_edges.csv", EDGE_FIELDS, result["edges"])
    _csv(directory / "country_network_branches.csv", BRANCH_FIELDS, result["branches"])
    _csv(
        directory / "country_network_node_uncertainty.csv",
        [
            "node_id",
            "country",
            "marginal_probability",
            "maximum_candidate",
            "annotated_representative",
            "is_tip",
        ],
        result["node_assignments"],
    )
    (directory / "country_network.json").write_text(json.dumps(result, indent=2) + "\n")
    _figure(directory, result)
    return result


def marginal_network_report_html(directory: Path) -> str:
    """A short report fragment with explicit reconstruction and uncertainty scope."""
    path = Path(directory) / "country_network.json"
    if not path.is_file():
        return "<p>Country network unavailable: no reconstruction summary was generated.</p>"
    result = json.loads(path.read_text())
    if result["status"] == "unavailable":
        return "<p>Country network unavailable: " + escape(result["reason"]) + "</p>"
    threshold = f"{100 * result['confidence_threshold']:g}%"
    rows = "".join(
        f"<tr><td>{escape(edge['source'])}</td><td>{escape(edge['target'])}</td>"
        f"<td>{edge['branches']}</td><td>{edge['high_endpoint_branches']}</td>"
        f"<td>{edge['uncertain_branches']}</td></tr>"
        for edge in result["edges"]
    )
    exclusions = len(result["excluded_tips"])
    exclusion_note = (
        f" {result['unknown_tips']} genomes have unknown location; {exclusions} excluded genomes are listed in the audit."
        if result["unknown_tips"] or exclusions
        else ""
    )
    return (
        f"<p>{result['changing_branches']} reconstructed location changes; "
        f"{result['uncertain_changing_branches']} involve uncertain location assignments.{exclusion_note}</p>"
        '<figure><a href="country_network.svg"><img src="country_network.svg" alt="Country network showing reconstructed '
        'country changes and uncertain endpoint assignments" style="max-width:100%;height:auto"></a>'
        "<figcaption>Circle counts show sampled genomes. Dashed arrows have uncertain location assignments; "
        f"solid arrows have both endpoint assignments ≥{threshold} with no tied maximum. "
        '<a href="country_network.svg">Full-resolution SVG</a> · <a href="country_network.png">PNG</a>.</figcaption></figure>'
        "<p>The chosen root and sampled genomes can change the arrows. Even a solid link is not evidence of direct transmission.</p>"
        "<details><summary>Country network data and method</summary>"
        f"<p>{escape(METHOD)}</p><p>{escape(CAVEAT)}</p><table><thead><tr><th>From</th><th>To</th><th>Branches</th>"
        "<th>Both endpoints meet threshold</th><th>Uncertain</th></tr></thead><tbody>"
        + rows
        + '</tbody></table><p><a href="country_network_nodes.csv">Country counts</a> · '
        '<a href="country_network_edges.csv">Arrow counts</a> · '
        '<a href="country_network_branches.csv">Every branch and endpoint confidence</a> · '
        '<a href="country_network_node_uncertainty.csv">All node state probabilities</a> · '
        '<a href="country_network.json">Full audit</a> · '
        '<a href="country_network.svg">SVG figure</a> · '
        '<a href="country_network.png">PNG figure</a></p></details>'
    )
