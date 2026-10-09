"""Offline country composition of frozen, annotated context catalogues.

No dates, tree tips, or live services participate in the public denominator.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import html
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

DEPTHS = (5, 6, 7)
CAUTION = (
    "Geographical composition of available sequenced records, affected by surveillance, "
    "submission coverage and selection. These figures do not estimate country prevalence, "
    "incidence, migration direction or transmission. Country metadata is not country of acquisition."
)
FIELDS = (
    "source",
    "snapshot",
    "scope",
    "filters",
    "scheme",
    "scheme_version",
    "prefix_depth",
    "prefix_key",
    "group_label",
    "assignment_category",
    "assignment_status",
    "cohort",
    "country",
    "count",
    "denominator",
    "known_country_n",
    "unknown_n",
    "percentage",
    "counting_unit",
    "raw_record_count",
)


def country_colour(country: str) -> str:
    """Stable across runs, cohorts and depths; Unknown and Other remain distinct."""
    if country == "Unknown":
        return "#9b9b9b"
    if country == "Other":
        return "#303030"
    import colorsys

    digest = hashlib.sha256(country.encode()).digest()
    hue = int.from_bytes(digest[:4], "big") / 2**32
    rgb = colorsys.hsv_to_rgb(hue, 0.48 + digest[4] / 255 * 0.22, 0.65 + digest[5] / 255 * 0.2)
    return "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)


def _source_id(row: Mapping[str, Any]) -> str:
    return str(row.get("source_genome_id") or row.get("genome_id") or row.get("id") or "")


def _aliases(row: Mapping[str, Any]) -> set[str]:
    aliases = set()
    for field, kind in (
        ("biosample", "biosample"),
        ("biosample_accession", "biosample"),
        ("assembly_accession", "assembly"),
        ("assembly", "assembly"),
    ):
        value = row.get(field)
        if value:
            aliases.add(f"{kind}:{str(value).strip().upper()}")
    for field, kind in (("assembly_accessions", "assembly"), ("run_accessions", "run")):
        for value in row.get(field, []) or []:
            aliases.add(f"{kind}:{str(value).strip().upper()}")
    for value in row.get("aliases", []) or []:
        aliases.add("accession:" + str(value).strip().upper())
    if row.get("sample_unit_id"):
        aliases.add("unit:" + str(row["sample_unit_id"]))
    if _source_id(row):
        aliases.add("source:" + str(row.get("source", "pathogenwatch")) + ":" + _source_id(row))
    return aliases


def _eligible(row: Mapping[str, Any]) -> bool:
    if "qc_pass" in row:
        return row["qc_pass"] is True
    if "qc_passed" in row:
        return row["qc_passed"] is True
    return str(row.get("qc_status", "")).lower() in {"pass", "passed"}


def _deduplicate(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Connected accession aliases handle assembly chains without counting them twice."""
    parent = list(range(len(rows)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for i, row in enumerate(rows):
        for alias in _aliases(row):
            if alias in seen:
                parent[find(i)] = find(seen[alias])
            seen[alias] = i
    clusters = defaultdict(list)
    for i, row in enumerate(rows):
        clusters[find(i)].append(row)
    units, audit = [], []
    for members in clusters.values():
        members.sort(key=lambda r: (str(r.get("source", "pathogenwatch")), _source_id(r)))
        unit = dict(members[0])
        aliases = sorted(set().union(*(_aliases(r) for r in members)))
        biosamples = [a for a in aliases if a.startswith("biosample:")]
        unit["_unit"] = (biosamples or aliases or ["unlinked:" + str(len(units))])[0]
        unit["_aliases"] = aliases
        unit["_source_ids"] = sorted(
            {
                identifier
                for r in members
                for identifier in (r.get("source_genome_ids") or [_source_id(r)])
            }
        )
        unit["_raw_n"] = sum(int(r.get("raw_genome_count", 1)) for r in members)
        countries = {str(r.get("country") or "Unknown") for r in members}
        if len(countries) > 1:
            unit["country"] = "Unknown"
        # Conflicting biological-sample assignments remain unresolved at every depth.
        signatures = {
            (
                str(r.get("cglin_raw", "")).replace(",", ".").replace("_", "."),
                str(r.get("cgst", "")),
                str(r.get("cglin_provisional", "")),
                str(r.get("cglin_scheme", "")),
                str(r.get("cglin_scheme_version", "")),
            )
            for r in members
        }
        assignment_conflict = len(signatures) > 1
        if assignment_conflict:
            unit["cglin_status"] = "conflict"
        for depth in DEPTHS:
            keys = {r.get(f"cglin_group_{depth}") for r in members if r.get(f"cglin_group_{depth}")}
            if assignment_conflict or len(keys) > 1:
                unit[f"cglin_group_{depth}"] = ""
                unit[f"cglin_status_{depth}"] = "conflict"
        for member in members:
            audit.append(
                {
                    "sample_unit": unit["_unit"],
                    "source_genome_id": _source_id(member),
                    "representative_source_id": _source_id(unit),
                    "sample_id": member.get("sample_id", ""),
                    "raw_country": member.get("country_raw", member.get("country", "")),
                    "normalised_country": member.get("country") or "Unknown",
                    "country_conflict": len(countries) > 1,
                    "cglin_conflict": assignment_conflict,
                    "biosample_conflict": len(biosamples) > 1,
                    "identity_resolved": bool(biosamples) or bool(unit.get("identity_resolved")),
                    "raw_records_in_unit": unit["_raw_n"],
                }
            )
        units.append(unit)
    return sorted(units, key=lambda r: r["_unit"]), audit


def geography_tables(
    catalogue_rows: Iterable[Mapping[str, Any]],
    *,
    selected_source_ids: Iterable[str] = (),
    focal_rows: Iterable[Mapping[str, Any]] = (),
    depths: Iterable[int] = DEPTHS,
    scope: Mapping[str, Any] | None = None,
) -> dict:
    """Return auditable full catalogue, selected and separate focal composition tables."""
    scope = dict(scope or {})
    rows = [dict(r) for r in catalogue_rows]
    eligible = [r for r in rows if _eligible(r)]
    units, audit = _deduplicate(eligible)
    selected = set(map(str, selected_source_ids))
    selected_units = [r for r in units if selected.intersection(r["_source_ids"])]
    focal_units, focal_audit = _deduplicate([dict(r) for r in focal_rows])
    public_aliases = set().union(*(set(r["_aliases"]) for r in units)) if units else set()
    overlap = []
    for unit in focal_units:
        matches = sorted(set(unit["_aliases"]) & public_aliases)
        if matches:
            overlap.append({"focal_sample_unit": unit["_unit"], "matched_aliases": matches})
    tables = []
    summaries = []
    for depth in depths:
        if depth not in DEPTHS:
            raise ValueError("Supported cgLIN prefix depths are 5, 6 and 7")
        for cohort, cohort_units in (
            ("public_catalogue", units),
            ("selected_context", selected_units),
            ("focal_survey", focal_units),
        ):
            groups = defaultdict(list)
            for unit in cohort_units:
                key = unit.get(f"cglin_group_{depth}") or ""
                status = unit.get(f"cglin_status_{depth}", unit.get("cglin_status", "missing"))
                valid = bool(key) and status in {"resolved", "provisional", "complete", "partial"}
                groups[
                    (
                        key if valid else "",
                        "lineage" if valid else "assignment_coverage",
                        "" if valid else str(status),
                        str(unit.get("cglin_scheme", "unknown")),
                        str(unit.get("cglin_scheme_version", "unknown")),
                        str(unit.get("source", "pathogenwatch")),
                    )
                ].append(unit)
            assigned = 0
            for (key, category, missing_status, scheme, version, source), members in sorted(
                groups.items()
            ):
                if category == "lineage":
                    assigned += len(members)
                counts = Counter(str(r.get("country") or "Unknown") for r in members)
                raw_counts = Counter()
                for member in members:
                    raw_counts[str(member.get("country") or "Unknown")] += member["_raw_n"]
                representative = members[0]
                label = key
                if key:
                    try:
                        parts = json.loads(key)
                        label = f"{parts[0]} / {parts[1]}: " + ".".join(map(str, parts[2]))
                    except (ValueError, TypeError, IndexError):
                        label = key
                else:
                    label = f"Unresolved: {missing_status} ({scheme} / {version}; {source})"
                for country, count in sorted(counts.items()):
                    tables.append(
                        dict(
                            zip(
                                FIELDS,
                                (
                                    source,
                                    scope.get(
                                        "snapshot", representative.get("snapshot", "unspecified")
                                    ),
                                    scope.get("description", "Frozen public same-ST catalogue"),
                                    scope.get(
                                        "filters",
                                        "Explicit QC pass; deduplicated; includes undated records",
                                    ),
                                    scheme,
                                    version,
                                    depth,
                                    key,
                                    label,
                                    category,
                                    "|".join(
                                        sorted(
                                            {
                                                str(
                                                    r.get(
                                                        f"cglin_status_{depth}",
                                                        r.get("cglin_status", "missing"),
                                                    )
                                                )
                                                for r in members
                                            }
                                        )
                                    ),
                                    cohort,
                                    country,
                                    count,
                                    len(members),
                                    len(members) - counts["Unknown"],
                                    counts["Unknown"],
                                    count / len(members) * 100,
                                    "deduplicated biological sample where accessions support identity; otherwise record unit",
                                    raw_counts[country],
                                ),
                            )
                        )
                    )
            summaries.append(
                {
                    "depth": depth,
                    "cohort": cohort,
                    "total_units": len(cohort_units),
                    "assigned_units": assigned,
                    "unresolved_units": len(cohort_units) - assigned,
                    "lineage_groups": sum(k[1] == "lineage" for k in groups),
                    "singleton_groups": sum(
                        k[1] == "lineage" and len(v) == 1 for k, v in groups.items()
                    ),
                }
            )
    return {
        "rows": tables,
        "summaries": summaries,
        "duplicate_audit": audit,
        "focal_duplicate_audit": focal_audit,
        "focal_overlap": overlap,
        "raw_catalogue_records": sum(int(r.get("raw_genome_count", 1)) for r in rows),
        "qc_eligible_records": sum(int(r.get("raw_genome_count", 1)) for r in eligible),
        "sample_units": len(units),
        "qc_excluded_records": sum(
            int(r.get("raw_genome_count", 1)) for r in rows if not _eligible(r)
        ),
        "identity_unresolved_units": sum(
            not (
                r.get("identity_resolved") or any(a.startswith("biosample:") for a in r["_aliases"])
            )
            for r in units
        ),
        "scope": scope,
    }


def _write_table(path: Path, rows: list[dict], fields: Iterable[str], delimiter=","):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def focused_geography_tables(result: Mapping[str, Any]) -> dict:
    """A bounded report view; full catalogue rows and denominators stay unchanged.

    Use the deepest level represented by a resolved focal assignment. Other
    focal assignments remain visibly unresolved at this level, not silently
    compared using a coarser prefix in the same composition denominator.
    """
    rows = result["rows"]
    available_depths = sorted({int(r["prefix_depth"]) for r in rows}, reverse=True)
    depth = next(
        (
            d
            for d in available_depths
            if any(
                r["cohort"] == "focal_survey"
                and r["prefix_depth"] == d
                and r["assignment_category"] == "lineage"
                for r in rows
            )
        ),
        None,
    )
    base_depth = min(available_depths) if available_depths else None
    counts = []
    names = {"focal_survey": "Focal samples", "selected_context": "Selected public context"}
    for cohort, label in names.items():
        country_counts = Counter()
        for row in rows:
            if row["cohort"] == cohort and row["prefix_depth"] == base_depth:
                country_counts[row["country"]] += row["count"]
        denominator = sum(country_counts.values())
        for country, count in sorted(country_counts.items()):
            counts.append(
                {
                    "cohort": cohort,
                    "cohort_label": label,
                    "country": country,
                    "count": count,
                    "denominator": denominator,
                    "percentage": count / denominator * 100 if denominator else 0,
                }
            )
    keys = sorted(
        {
            r["prefix_key"]
            for r in rows
            if r["cohort"] == "focal_survey"
            and r["prefix_depth"] == depth
            and r["assignment_category"] == "lineage"
        }
    )
    labels = {key: f"Group {index + 1}" for index, key in enumerate(keys)}
    lookup, public = [], []
    for key in keys:
        focal = [
            r
            for r in rows
            if r["cohort"] == "focal_survey"
            and r["prefix_depth"] == depth
            and r["prefix_key"] == key
        ]
        contexts = [
            r
            for r in rows
            if r["cohort"] == "selected_context"
            and r["prefix_depth"] == depth
            and r["prefix_key"] == key
        ]
        matching = [
            r
            for r in rows
            if r["cohort"] == "public_catalogue"
            and r["prefix_depth"] == depth
            and r["prefix_key"] == key
        ]
        representative = focal[0]
        lookup.append(
            {
                "group_id": labels[key],
                "prefix_depth": depth,
                "prefix_key": key,
                "full_prefix_label": representative["group_label"],
                "scheme": representative["scheme"],
                "scheme_version": representative["scheme_version"],
                "focal_n": sum(r["count"] for r in focal),
                "selected_context_n": sum(r["count"] for r in contexts),
                "public_n": sum(r["count"] for r in matching),
            }
        )
        for row in matching:
            public.append(dict(row, group_id=labels[key]))
    chosen_summaries = [dict(s) for s in result["summaries"] if s["depth"] == depth]
    fallback_summaries = [dict(s) for s in result["summaries"] if s["depth"] == base_depth]
    totals = {s["cohort"]: s["total_units"] for s in fallback_summaries}
    unresolved = [
        dict(r)
        for r in rows
        if r["prefix_depth"] == (depth or base_depth)
        and r["assignment_category"] == "assignment_coverage"
    ]
    return {
        "depth": depth,
        "country_counts": counts,
        "public_groups": public,
        "group_lookup": lookup,
        "assignment_coverage": unresolved,
        "audit": {
            "scope": dict(result.get("scope", {})),
            "chosen_depth": depth,
            "public_catalogue_n": result["sample_units"],
            "focal_n": totals.get("focal_survey", 0),
            "selected_context_n": totals.get("selected_context", 0),
            "tree_participants_n": totals.get("focal_survey", 0)
            + totals.get("selected_context", 0),
            "focal_public_overlap_n": len(result["focal_overlap"]),
            "matched_public_n": sum(g["public_n"] for g in lookup),
            "groups_without_public_records": [g["group_id"] for g in lookup if not g["public_n"]],
            "assignment_summaries": chosen_summaries or fallback_summaries,
            "counting_unit": "deduplicated sample units; not patients or infections",
            "denominator_rule": "Each public group has its own denominator including Unknown; focal and selected counts are separate.",
            "depth_rule": "Deepest resolved focal prefix, 7 then 6 then 5; unresolved focal assignments remain outside these groups.",
            "interpretation": "Shared cgLIN group membership is not evidence of exact nearest relatives or transmission.",
        },
    }


def _focused_plot(rows: list[dict], path: Path, *, proportions: bool):
    """Named country bars with readable n/% labels, kept out of the HTML DOM."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    with plt.rc_context({"font.family": "DejaVu Serif", "svg.fonttype": "none", "font.size": 14}):
        if proportions:
            groups = list(dict.fromkeys(r["group_id"] for r in rows))
            panels = [[r for r in rows if r["group_id"] == group] for group in groups]
        else:
            panels = [rows]
        countries_per_panel = [len({r["country"] for r in panel}) for panel in panels]
        fig, axes = plt.subplots(
            len(panels),
            1,
            squeeze=False,
            figsize=(8, max(4.6, sum(countries_per_panel) * 0.48 + len(panels) * 2)),
        )
        for ax, panel in zip(axes.flat, panels):
            countries = sorted(
                {r["country"] for r in panel},
                key=lambda country: (
                    country == "Unknown",
                    -sum(r["count"] for r in panel if r["country"] == country),
                    country,
                ),
            )
            if proportions:
                by_country = {r["country"]: r for r in panel}
                values = [by_country[c]["percentage"] for c in countries]
                bars = ax.barh(
                    range(len(countries)),
                    values,
                    color=[country_colour(c) for c in countries],
                    height=0.65,
                )
                for bar, country in zip(bars, countries):
                    row = by_country[country]
                    ax.text(
                        bar.get_width() + 1,
                        bar.get_y() + bar.get_height() / 2,
                        f"{row['count']:,} ({row['percentage']:.1f}%)",
                        va="center",
                        fontsize=14,
                    )
                ax.set_xlim(0, max(25, max(values, default=0) * 1.55))
                ax.set_xlabel("Share of this public group (%)\nUnknown included in denominator")
                first = panel[0]
                ax.set_title(
                    f"{first['group_id']} · prefix depth {first['prefix_depth']} · public N={first['denominator']:,}\n"
                    f"Known country={first['known_country_n']:,} · Unknown={first['unknown_n']:,}",
                    loc="left",
                    fontsize=14,
                )
            else:
                offsets = {"focal_survey": -0.18, "selected_context": 0.18}
                colours = {"focal_survey": "#303a46", "selected_context": "#b78151"}
                for cohort, offset in offsets.items():
                    cohort_rows = [r for r in panel if r["cohort"] == cohort]
                    by_country = {r["country"]: r["count"] for r in cohort_rows}
                    denominator = sum(by_country.values())
                    label = (
                        "Focal samples" if cohort == "focal_survey" else "Selected public context"
                    )
                    values = [by_country.get(c, 0) for c in countries]
                    bars = ax.barh(
                        [i + offset for i in range(len(countries))],
                        values,
                        height=0.32,
                        color=colours[cohort],
                        label=f"{label} (N={denominator:,})",
                    )
                    for bar, value in zip(bars, values):
                        if value:
                            ax.text(
                                value + 0.04,
                                bar.get_y() + bar.get_height() / 2,
                                str(value),
                                va="center",
                                fontsize=14,
                            )
                ax.set_xlim(0, max(1, max((r["count"] for r in panel), default=0)) * 1.25)
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax.set_xlabel("Sample units (n)")
                ax.set_title("Countries represented in the tree cohorts", loc="left", fontsize=14)
                ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -0.28), fontsize=14)
            ax.set_yticks(range(len(countries)), countries)
            ax.invert_yaxis()
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="x", color="#dddddd", linewidth=0.5, zorder=0)
            ax.set_axisbelow(True)
        fig.tight_layout(pad=1.6)
        fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
        fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight")
        plt.close(fig)


def generate_focused_geography(result: Mapping[str, Any], output_dir: str | Path) -> dict:
    """Save a lightweight main-report summary plus separate source tables/audit."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    focused = focused_geography_tables(result)
    for previous in output_dir.glob("focused_*.svg"):
        previous.unlink()
    for previous in output_dir.glob("focused_*.png"):
        previous.unlink()
    files = []
    for name, rows, fields in (
        (
            "focused_country_counts.csv",
            focused["country_counts"],
            ("cohort", "cohort_label", "country", "count", "denominator", "percentage"),
        ),
        ("focused_public_groups.csv", focused["public_groups"], (*FIELDS, "group_id")),
        (
            "focused_group_lookup.csv",
            focused["group_lookup"],
            (
                "group_id",
                "prefix_depth",
                "prefix_key",
                "full_prefix_label",
                "scheme",
                "scheme_version",
                "focal_n",
                "selected_context_n",
                "public_n",
            ),
        ),
        ("focused_assignment_coverage.csv", focused["assignment_coverage"], FIELDS),
    ):
        _write_table(output_dir / name, rows, fields)
        files.append(str(output_dir / name))
    audit_path = output_dir / "focused_scope_audit.json"
    audit_path.write_text(json.dumps(focused["audit"], indent=2) + "\n")
    files.append(str(audit_path))
    audit = focused["audit"]
    fragment = [
        '<section id="context-geography">',
        f"<p>The tree includes {audit['focal_n']:,} focal samples and {audit['selected_context_n']:,} selected public comparisons "
        f"({audit['tree_participants_n']:,} genomes in total). The wider catalogue contains {audit['public_catalogue_n']:,} "
        "quality-checked, deduplicated public sample units; these are a separate population, not extra tree participants.</p>",
    ]
    scope = result.get("scope", {})
    for name, rows, proportions, caption in (
        (
            "focused_country_counts",
            focused["country_counts"],
            False,
            "Figure 1. Countries of the focal samples and selected public context, with separate cohort counts. All named countries and Unknown are retained.",
        ),
        (
            "focused_public_groups",
            focused["public_groups"],
            True,
            "Figure 2. Country proportions among available public records within each focal-matching cgLIN group. Labels show count and percentage; Unknown remains in each group denominator. Focal counts are not added to public proportions.",
        ),
    ):
        if rows:
            if proportions:
                fragment.append(
                    "<p>Public proportions describe genetic groups matching your focal samples; "
                    "they do not identify exact nearest relatives.</p>"
                )
            _focused_plot(rows, output_dir / name, proportions=proportions)
            files.extend(str(output_dir / (name + suffix)) for suffix in (".svg", ".png"))
            fragment.append(
                f'<figure><a href="{name}.svg" target="_blank"><img src="{name}.svg" '
                f'alt="{html.escape(caption)}" loading="lazy" style="width:100%;height:auto"></a>'
                f'<figcaption>{html.escape(caption)} <a href="{name}.svg">Full-resolution SVG</a> · '
                f'<a href="{name}.png">PNG</a></figcaption></figure>'
            )
    if scope:
        fragment.append(
            "<details><summary>Catalogue source, snapshot and quality filters</summary><p>Available public records: "
            + html.escape(str(scope.get("description", "Frozen same-ST catalogue")))
            + ". Snapshot: "
            + html.escape(str(scope.get("snapshot", "unspecified")))
            + ". "
            + html.escape(
                str(scope.get("filters", "QC passing; deduplicated; undated records retained"))
            )
            + ".</p></details>"
        )
    if focused["depth"] is None:
        fragment.append(
            "<p>Public group composition is unavailable: no comparable focal cgLIN assignment resolves at depths 7, 6 or 5. "
            "Missing and unresolved assignments remain in the coverage table; no lineage group has been inferred.</p>"
        )
    else:
        focal_summary = next(
            s for s in audit["assignment_summaries"] if s["cohort"] == "focal_survey"
        )
        fragment.append(
            f"<p>Comparison uses full-prefix depth {focused['depth']}; {focal_summary['assigned_units']:,} focal sample units resolve "
            f"and {focal_summary['unresolved_units']:,} remain unresolved at this depth. Prefix depths are not SNP cutoffs.</p>"
        )
        if audit["groups_without_public_records"]:
            fragment.append(
                "<p>No eligible public records were available for "
                + html.escape(", ".join(audit["groups_without_public_records"]))
                + ".</p>"
            )
    if focused["group_lookup"]:
        fragment.append(
            "<details><summary>Group identifiers, schemes and separate cohort totals</summary>"
            "<table><thead><tr><th>Group</th><th>Full prefix, scheme and version</th><th>Public N</th><th>Focal N</th><th>Selected N</th></tr></thead><tbody>"
        )
        for group in focused["group_lookup"]:
            fragment.append(
                "<tr><td>"
                + html.escape(group["group_id"])
                + "</td><td>"
                + html.escape(group["full_prefix_label"])
                + f"</td><td>{group['public_n']:,}</td><td>{group['focal_n']:,}</td><td>{group['selected_context_n']:,}</td></tr>"
            )
        fragment.append(
            "</tbody></table><p>The scheme version is retained exactly, including Unknown where unavailable. "
            "A Pathogenwatch source ID identifies a source genome record; it is not a lineage or patient identifier.</p></details>"
        )
    fragment.append(
        f"<p>Focal sample units also present in the public catalogue: {audit['focal_public_overlap_n']:,}. "
        "Their focal counts remain separate and do not increase public proportions.</p>"
    )
    fragment.append("<p>" + html.escape(CAUTION) + "</p>")
    fragment.append(
        '<p><a href="focused_country_counts.csv">Cohort country counts</a> · '
        '<a href="focused_public_groups.csv">Focal-matching public proportions</a> · '
        '<a href="focused_group_lookup.csv">Full group lookup</a> · '
        '<a href="focused_assignment_coverage.csv">Unresolved assignment coverage</a> · '
        '<a href="focused_scope_audit.json">Scope audit</a></p>'
    )
    fragment.append(
        "<details><summary>Explore the complete country atlas</summary><p>"
        '<a href="index.html">Open all catalogue groups, depths and cohort figures</a> · '
        '<a href="country_composition.csv">Download the complete audited country table</a>.'
        "</p></details></section>"
    )
    (output_dir / "focused_fragment.html").write_text("\n".join(fragment))
    files.append(str(output_dir / "focused_fragment.html"))
    return {"outputs": files, "focused": focused, "focused_report_html": "\n".join(fragment)}


def _plot(rows: list[dict], path: Path, title: str, percentage: bool, caption: str):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    countries = sorted({r["country"] for r in rows})
    if len(countries) > 13:
        raise ValueError("Plot country grouping was not applied")
    groups = list(
        dict.fromkeys(
            (
                r["prefix_key"],
                r["group_label"],
                r["denominator"],
                r["known_country_n"],
                r["unknown_n"],
            )
            for r in rows
        )
    )
    fig, ax = plt.subplots(figsize=(14, max(5, len(groups) * 0.48 + 3)))
    lookup = {(r["prefix_key"], r["group_label"], r["country"]): r for r in rows}
    left = [0.0] * len(groups)
    for country in countries:
        values = [
            lookup.get((key, label, country), {}).get("percentage" if percentage else "count", 0)
            for key, label, *_ in groups
        ]
        ax.barh(
            range(len(groups)),
            values,
            left=left,
            label=country,
            color=country_colour(country),
            edgecolor="white",
            linewidth=0.4,
        )
        left = [a + b for a, b in zip(left, values)]
    ax.set_yticks(
        range(len(groups)),
        [
            f"{label}  N={n}; known={known}; Unknown={unknown}"
            + (" (singleton)" if n == 1 else " (small)" if n < 5 else "")
            for _, label, n, known, unknown in groups
        ],
        fontsize=9,
    )
    ax.invert_yaxis()
    ax.set_xlabel(
        "Within-group percentage: n/N, including Unknown (%)"
        if percentage
        else "Deduplicated sample units (n); Unknown included in N"
    )
    if percentage:
        ax.set_xlim(0, 100)
    ax.set_title(title, fontsize=13, pad=15)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=min(5, len(countries)),
        frameon=False,
        fontsize=9,
    )
    fig.text(0.02, 0.015, caption, fontsize=9, va="bottom", wrap=True)
    footer = max(0.28, min(0.44, 2.2 / fig.get_size_inches()[1]))
    fig.subplots_adjust(left=0.35, right=0.98, bottom=footer, top=0.88)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def generate_context_geography(
    catalogue_rows: Iterable[Mapping[str, Any]],
    output_dir: str | Path,
    *,
    selected_source_ids: Iterable[str] = (),
    focal_rows: Iterable[Mapping[str, Any]] = (),
    depths: Iterable[int] = DEPTHS,
    scope: Mapping[str, Any] | None = None,
) -> dict:
    """Write CSV/TSV, JSON audit, figures, and self-contained embeddable HTML fragment."""
    catalogue_rows = list(catalogue_rows)
    selected_source_ids = list(selected_source_ids)
    focal_rows = list(focal_rows)
    depths = tuple(depths)
    result = geography_tables(
        catalogue_rows,
        selected_source_ids=selected_source_ids,
        focal_rows=focal_rows,
        depths=depths,
        scope=scope,
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Reusing an output directory must not archive plots from a previous cohort/snapshot.
    for cohort in ("public_catalogue", "selected_context", "focal_survey"):
        for suffix in ("svg", "png"):
            for previous in output_dir.glob(f"{cohort}_depth_*.{suffix}"):
                previous.unlink()
    for name in ("country_composition.html", "fragment.html", "index.html"):
        (output_dir / name).unlink(missing_ok=True)
    _write_table(output_dir / "country_composition.csv", result["rows"], FIELDS)
    _write_table(output_dir / "country_composition.tsv", result["rows"], FIELDS, "\t")
    matrix = {}
    countries = sorted({r["country"] for r in result["rows"]})
    matrix_fields = [
        f for f in FIELDS if f not in {"country", "count", "percentage", "raw_record_count"}
    ]
    for item in result["rows"]:
        key = tuple(item[f] for f in matrix_fields)
        if key not in matrix:
            matrix[key] = {f: item[f] for f in matrix_fields}
            matrix[key].update({country: 0 for country in countries})
        matrix[key][item["country"]] = item["count"]
    _write_table(
        output_dir / "country_composition_matrix.csv",
        list(matrix.values()),
        matrix_fields + countries,
    )
    (output_dir / "country_composition_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    outputs = [
        str(output_dir / name)
        for name in (
            "country_composition.csv",
            "country_composition.tsv",
            "country_composition_matrix.csv",
            "country_composition_audit.json",
        )
    ]
    fragments = [
        '<section id="context-geography"><h2>Country composition within cgLIN groups</h2>',
        "<p>" + html.escape(CAUTION) + "</p>",
        "<p>Public denominators include QC-passing undated records before selection. "
        "Focal denominators are separate; public overlap is reported in the audit.</p>",
    ]
    if not result["rows"]:
        fragments.append("<p>No eligible annotated context or focal records are available.</p>")
    fragments.append(
        f"<p>Focal sample units overlapping public catalogue: {len(result['focal_overlap'])}; "
        "their survey counts remain separate from public proportions.</p>"
    )
    for name in (
        "country_composition.csv",
        "country_composition.tsv",
        "country_composition_matrix.csv",
        "country_composition_audit.json",
    ):
        encoded = base64.b64encode((output_dir / name).read_bytes()).decode()
        fragments.append(
            f'<p><a download="{name}" href="data:application/octet-stream;base64,{encoded}">Download {name}</a></p>'
        )
    country_totals = Counter()
    for item in result["rows"]:
        if item["cohort"] == "public_catalogue" and item["prefix_depth"] == min(depths):
            country_totals[item["country"]] += item["count"]
    if not country_totals:
        for item in result["rows"]:
            if item["prefix_depth"] == min(depths):
                country_totals[item["country"]] += item["count"]
    visible_countries = {
        c
        for c, _ in sorted(country_totals.items(), key=lambda pair: (-pair[1], pair[0]))
        if c != "Unknown"
    }
    visible_countries = set(sorted(visible_countries, key=lambda c: (-country_totals[c], c))[:11])
    visible_countries.add("Unknown")
    for summary in result["summaries"]:
        cohort, depth = summary["cohort"], summary["depth"]
        cohort_countries = {r["country"] for r in result["rows"] if r["cohort"] == cohort}
        show_all = cohort != "public_catalogue" and len(cohort_countries) <= 12
        panel_countries = cohort_countries if show_all else visible_countries
        display_rule = (
            "All named countries shown for this bounded cohort."
            if show_all
            else "Top 11 named countries by public-catalogue count + Unknown; remaining countries = Other. Full named-country table downloadable."
        )
        denominator_note = {
            "public_catalogue": "Full QC-passing deduplicated public catalogue; undated records retained.",
            "selected_context": "Selected contexts after date eligibility, bounded pool, downloads and SKA selection.",
            "focal_survey": "Separate focal sample denominator; public overlaps are not added to public proportions.",
        }[cohort]
        if depth != 5:
            fragments.append(
                "<details><summary>"
                + html.escape(
                    f"{cohort}: explore prefix depth {depth}; N={summary['total_units']}, "
                    f"unresolved={summary['unresolved_units']}"
                )
                + "</summary>"
            )
        fragments.append("<h3>" + html.escape(f"{cohort}: full-prefix depth {depth}") + "</h3>")
        fragments.append(
            f"<p>Total N={summary['total_units']}; resolved N={summary['assigned_units']}; "
            f"unresolved N={summary['unresolved_units']}. Prefix depths are not SNP cutoffs.</p>"
        )
        for category in ("lineage", "assignment_coverage"):
            panel = [
                r
                for r in result["rows"]
                if r["cohort"] == cohort
                and r["prefix_depth"] == depth
                and r["assignment_category"] == category
            ]
            if not panel:
                if category == "lineage":
                    fragments.append("<p>No resolved comparable cgLIN groups at this depth.</p>")
                continue
            # Fixed page size retains singletons without making labels unreadable.
            panel.sort(
                key=lambda r: (-r["denominator"], r["prefix_key"], r["group_label"], r["country"])
            )
            labels = list(dict.fromkeys((r["prefix_key"], r["group_label"]) for r in panel))
            for page in range(0, len(labels), 24):
                page_labels = set(labels[page : page + 24])
                page_rows = [r for r in panel if (r["prefix_key"], r["group_label"]) in page_labels]
                for mode in ("counts", "percentages"):
                    stem = f"{cohort}_depth_{depth}_{category}_{mode}_{page // 24 + 1}"
                    visual = {}
                    for item in page_rows:
                        item = dict(item)
                        if item["country"] not in panel_countries:
                            item["country"] = "Other"
                        merge_key = (item["prefix_key"], item["group_label"], item["country"])
                        if merge_key in visual:
                            for field in ("count", "percentage", "raw_record_count"):
                                visual[merge_key][field] += item[field]
                        else:
                            visual[merge_key] = item
                    path = output_dir / stem
                    title = f"{cohort} — depth {depth} — {category.replace('_', ' ')} — {mode}"
                    caption = (
                        "Scope: "
                        + str(result["scope"].get("description", "Frozen public same-ST catalogue"))
                        + "\nSnapshot: "
                        + str(result["scope"].get("snapshot", "unspecified"))
                        + ". Source: "
                        + ", ".join(sorted({str(r["source"]) for r in page_rows}))
                        + "\nFilters: "
                        + str(
                            result["scope"].get(
                                "filters",
                                "Explicit QC pass; deduplicated; undated records included",
                            )
                        )
                        + "\nDenominator: "
                        + denominator_note
                        + "\nDisplay: "
                        + display_rule
                        + "\n"
                        + CAUTION
                    )
                    _plot(list(visual.values()), path, title, mode == "percentages", caption)
                    outputs.extend([str(path.with_suffix(".svg")), str(path.with_suffix(".png"))])
                    svg = path.with_suffix(".svg").read_text()
                    fragments.append("<svg" + svg.split("<svg", 1)[1])
        if depth != 5:
            fragments.append("</details>")
    fragments.append("</section>")
    result["outputs"] = outputs
    result["report_html"] = "\n".join(fragments)
    for name in ("country_composition.html", "fragment.html"):
        (output_dir / name).write_text(result["report_html"])
        result["outputs"].append(str(output_dir / name))
    document = (
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        "<title>cgLIN country composition</title><style>"
        "body{font:16px system-ui,sans-serif;margin:2rem;max-width:1500px}"
        "svg{width:100%;height:auto;display:block;margin:1rem 0}"
        "</style><body>" + result["report_html"] + "</body></html>"
    )
    (output_dir / "index.html").write_text(document)
    result["outputs"].append(str(output_dir / "index.html"))
    focused_result = (
        result
        if set(depths) == set(DEPTHS)
        else geography_tables(
            catalogue_rows,
            selected_source_ids=selected_source_ids,
            focal_rows=focal_rows,
            depths=DEPTHS,
            scope=scope,
        )
    )
    focused = generate_focused_geography(focused_result, output_dir)
    result["outputs"].extend(focused["outputs"])
    result["focused"] = focused["focused"]
    result["focused_report_html"] = focused["focused_report_html"]
    return result
