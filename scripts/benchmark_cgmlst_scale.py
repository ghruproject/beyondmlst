"""Run a labelled synthetic full-profile distance/NJ/PCoA/network resource benchmark."""

import argparse
import json
import resource
import time
import sys
import platform
from pathlib import Path
import numpy as np
from chronoclade.cgmlst.scale import analyse_profiles_scale

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--samples", type=int, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--countries", type=int, default=5, help="Reported categories including one unknown"
)
parser.add_argument(
    "--country-independent",
    action="store_true",
    help="Decouple reported country from genetic template",
)
args = parser.parse_args()
if args.samples < 20 or args.countries < 2:
    parser.error("Need at least 20 samples and two reported country categories")
size = args.samples
output = args.output
rng = np.random.default_rng(173)
countries = (
    ["UK", "India", "France", "Spain"]
    if args.countries == 5
    else [f"Synthetic country {i}" for i in range(args.countries - 1)]
)
countries.append("Unknown")
loci = [f"locus_{i:04d}" for i in range(629)]
rows = []
for i in range(size):
    group = i % 20
    calls = np.full(629, group + 1)
    changed = rng.choice(629, 3, replace=False)
    calls[changed] = rng.integers(21, 201, size=3)
    rows.append(
        dict(
            sample_id=f"SYNTHETIC_{i:05d}",
            role="input" if i < 10 else "context",
            origin="local" if i < 10 else "context",
            species="Synthetic categorical profiles",
            country=countries[(i // 20 if args.country_independent else i) % len(countries)],
            cgmlst_scheme="synthetic629",
            cgmlst_scheme_version="simulation-1",
            cgmlst_loci=loci,
            cgmlst_profile=dict(zip(loci, calls.tolist())),
            collection_date=str(2018 + i % 8),
        )
    )
started = time.perf_counter()
result = analyse_profiles_scale(
    rows, output=output, bootstrap_replicates=0, rapidnj_memory_mb=2048, batch_size=128
)
wall = time.perf_counter() - started
usage = resource.getrusage(resource.RUSAGE_SELF)
children = resource.getrusage(resource.RUSAGE_CHILDREN)
record = dict(
    schema="chronoclade.cgmlst.synthetic_scale_benchmark",
    sample_count=size,
    locus_count=629,
    synthetic=True,
    seed=173,
    model=f"20 categorical templates; three independently changed loci/profile; {args.countries} reported countries including unknown; ten inputs; eight collection years; country independent of template={args.country_independent}",
    hardware=dict(platform=platform.platform(), machine=platform.machine()),
    wall_seconds=wall,
    peak_self_rss_bytes=usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024,
    peak_child_rss_bytes=children.ru_maxrss
    if sys.platform == "darwin"
    else children.ru_maxrss * 1024,
    backend=result["scale_backend"],
    cohorts=[
        dict(
            profile_count=len(c["sample_ids"]),
            seconds=c["scale_runtime_seconds"],
            nj_backend=c.get("nj_backend", {}),
            pcoa_solver=c["pcoa_solver"],
        )
        for c in result["cohorts"]
    ],
    network=[
        dict(
            profile_count=n["sample_count"],
            seconds=n["scale_runtime_seconds"],
            country_count=len(n["nodes"]),
            root_count=n["root_count"],
            uncertainty_scope=n["uncertainty_scope"],
        )
        for n in result["location_network"]
    ],
    comparable_pairs=result["coverage"]["comparable_pairs"],
    available_profiles=result["coverage"]["available_profiles"],
    ordination_csv_rows=sum(
        len(Path(c["pcoa_csv"]).read_text().splitlines()) - 1 for c in result["cohorts"]
    ),
    tree_tip_count=sum(
        len(__import__("Bio").Phylo.read(c["tree_path"], "newick").get_terminals())
        for c in result["cohorts"]
    ),
    output_disk_bytes=sum(p.stat().st_size for p in output.rglob("*") if p.is_file()),
    interpretation="Computational resource/correctness exercise on labelled synthetic inputs; not evidence of biological adequacy, real-world cohort comparability, or large-country network behaviour.",
)
(output / "benchmark.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
