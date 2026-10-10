"""Run one profile analysis per input ST/CG dataset and retain its full context."""

from collections import Counter, defaultdict
from copy import deepcopy
import csv
import json
from pathlib import Path

from chronoclade.adaptive_context import annotate_cglin_datasets
from chronoclade.errors import WorkflowError
from chronoclade.metadata import slugify_lineage


def route_datasets(inputs):
    """Use CG dataset names for downstream stages, preserving the original ST."""
    result = deepcopy(inputs)
    for name in ("queries", "context", "catalogue_rows"):
        for row in result.get(name, []):
            normalized = " ".join(row.get("species", "").replace("_", " ").split()).lower()
            row["species"] = normalized[:1].upper() + normalized[1:]
    # Public fallback rows retain their original ST pathway. Only explicit user
    # comparisons need CG assignment here; selected adaptive rows are labelled.
    for index, row in enumerate(result["context"]):
        if row.get("provided_context") and not row.get("analysis_dataset"):
            result["context"][index] = annotate_cglin_datasets([row], result["queries"])[0]
    result["catalogue_rows"] = annotate_cglin_datasets(
        result.get("catalogue_rows", []), result["queries"]
    )
    if result["provenance"].get("adaptive_context_selection", {}).get("datasets"):
        query_ids = {r.get("source_genome_id") for r in result["queries"] if r.get("source_genome_id")}
        result["catalogue_rows"] = [r for r in result["catalogue_rows"]
                                    if r.get("source_genome_id") not in query_ids]
    for name in ("queries", "context", "catalogue_rows"):
        for row in result.get(name, []):
            if row.get("analysis_dataset"):
                row["input_lineage"] = row["lineage"]
                row["lineage"] = row["analysis_dataset"]
    return result


def geography(rows):
    counts = Counter((r.get("origin", "context"), r.get("country") or "Unknown",
                      r.get("region") or "Unknown", r.get("nuts2") or "Unknown") for r in rows)
    return [dict(origin=k[0], country=k[1], region=k[2], nuts2=k[3], count=n)
            for k, n in sorted(counts.items())]


def run_fast_datasets(inputs, output, *, seed, bootstrap_replicates,
                      distance_threshold, tree_limit=80, context_size=50,
                      nearest_per_query=3, include_genomes=None):
    from chronoclade.profile_analysis import analyse_profiles, set_tree_display_samples
    from chronoclade.metadata_dates import date_interval as _date_interval
    from chronoclade.profile_report import write_profile_report, write_fast_group_index
    from chronoclade.report import write_supporting_bundle
    from chronoclade.selections import select_assembly_context

    groups = defaultdict(list)
    for row in inputs["queries"]:
        groups[(row["species"], row["lineage"])].append(row)
    selectable = [r for r in inputs["context"] if (r["species"], r["lineage"]) in groups]
    if len({r["sample_id"] for r in inputs["queries"] + selectable}) != len(inputs["queries"] + selectable):
        raise WorkflowError("Duplicate sample identifiers in query lineage pools")
    requested = {}
    for identifier in include_genomes or []:
        matches = [r for r in selectable if identifier in {
            r["sample_id"], str(r.get("source_genome_id", "")), str(r.get("accession", ""))}]
        if len(matches) != 1:
            raise WorkflowError(
                f"Requested context genome {identifier!r} is not uniquely in a query lineage pool"
            )
        requested[identifier] = (matches[0]["species"], matches[0]["lineage"])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    reports, results, pairs, selections = [], [], [], []
    for key, queries in sorted(groups.items()):
        context = [r for r in inputs["context"] if (r["species"], r["lineage"]) == key]
        catalogue = [r for r in inputs.get("catalogue_rows", [])
                     if (r["species"], r["lineage"]) == key]
        directory = output if len(groups) == 1 else output / slugify_lineage(*key)
        analysis = analyse_profiles(queries + context, output=directory, seed=seed,
                                    bootstrap_replicates=bootstrap_replicates,
                                    distance_threshold=distance_threshold,
                                    tree_limit=tree_limit)
        selected, selection = select_assembly_context(
            queries, context, analysis, size=context_size,
            nearest_per_query=nearest_per_query,
            include=[identifier for identifier, dataset in requested.items() if dataset == key],
            seed=seed,
        )
        selection.update(
            selected_context_ids=[r["sample_id"] for r in selected],
            selected_sample_ids=[r["sample_id"] for r in queries + selected],
            available_context_ids=sorted(r["sample_id"] for r in context),
            query_ids=[r["sample_id"] for r in queries],
        )
        set_tree_display_samples(analysis, queries + context, selection["selected_sample_ids"])
        analysis["shared_selection"] = selection
        selections.append(dict(species=key[0], lineage=key[1], **selection))
        selection_path = directory / "context_selection.json"
        selection_path.write_text(json.dumps(selection, indent=2) + "\n")
        analysis.setdefault("paths", {})["context_selection"] = str(selection_path)
        analysis.update(species=key[0], lineage=key[1],
                        metadata_geography=geography(queries + context),
                        public_catalogue_geography=geography(catalogue))
        provenance = deepcopy(inputs["provenance"])
        provenance["shared_selection"] = selection
        provenance["coverage"] = {
            label: {"total": len(rows),
                    "profiles_available": sum(bool(r.get("cgmlst_profile") or r.get("cgmlst_novel_alleles")) for r in rows)}
            for label, rows in (("queries", queries), ("context", context))
        }
        adaptive = provenance.get("adaptive_context_selection")
        if adaptive:
            adaptive["datasets"] = [d for d in adaptive["datasets"] if d["analysis_dataset"] == key[1]]
            for dataset in adaptive["datasets"]:
                provenance["eligible_public_context_count"] = dataset["public_cg_pool_count"]
                provenance["bounded_public_context_count"] = len(context)
                provenance["catalogue_metadata_count"] = len(catalogue)
                provenance.pop("profile_limit", None)
                provenance.pop("context_pool_selection", None)
                dataset["selected_context_geography"] = [dict(r, denominator=len(context)) for r in geography(context)]
                years = Counter(str((_date_interval(r) or {}).get("year", "Unknown")) for r in context)
                dataset["selected_context_years"] = [dict(year=y, count=n, denominator=len(context))
                                                     for y, n in sorted(years.items())]
            selection_path = directory / "adaptive_context_selection.json"
            selection_path.write_text(json.dumps(adaptive, indent=2) + "\n")
            analysis.setdefault("paths", {})["adaptive_context_selection"] = str(selection_path)
        (directory / "profile_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        (directory / "profile_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
        (directory / "metadata_geography.json").write_text(json.dumps({
            "analysis_records": analysis["metadata_geography"],
            "public_catalogue": analysis["public_catalogue_geography"]}, indent=2) + "\n")
        report = write_profile_report(analysis, directory=directory, provenance=provenance)
        write_supporting_bundle(directory)
        reports.append(dict(species=key[0], lineage=key[1], input_count=len(queries),
                            context_count=len(context), report=str(report)))
        results.append(analysis)
        pair_path = analysis.get("paths", {}).get("pairwise_distances")
        if pair_path and Path(pair_path).is_file():
            with Path(pair_path).open() as handle:
                pairs.extend(csv.DictReader(handle))
    if len(results) == 1:
        results[0]["context_selections"] = selections
        (output / "profile_analysis.json").write_text(json.dumps(results[0], indent=2) + "\n")
        write_supporting_bundle(output)
        return results[0], Path(reports[0]["report"]), reports
    # Keep complete-pool distance evidence together for downstream inspection;
    # selections and figures remain separate per dataset.
    combined = {"nearest_neighbours": [], "cohorts": [], "genetic_groups": [], "paths": {},
                "context_selections": selections}
    for result in results:
        combined["nearest_neighbours"].extend(result.get("nearest_neighbours", []))
        combined["cohorts"].extend(result.get("cohorts", []))
        combined["genetic_groups"].extend(result.get("genetic_groups", []))
    if pairs:
        pair_path = output / "pairwise_distances.csv"
        with pair_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(pairs[0]))
            writer.writeheader()
            writer.writerows(pairs)
        combined["paths"]["pairwise_distances"] = str(pair_path)
    selection_path = output / "context_selection.json"
    selection_path.write_text(json.dumps(selections, indent=2) + "\n")
    combined["paths"]["context_selection"] = str(selection_path)
    report = write_fast_group_index(output, reports,
                                    provenance=dict(inputs["provenance"], shared_selection=selections))
    (output / "profile_analysis.json").write_text(json.dumps(combined, indent=2) + "\n")
    return combined, report, reports
