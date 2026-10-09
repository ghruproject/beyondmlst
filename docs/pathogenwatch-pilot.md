# Reproducible ST147 pilot on a laptop and SLURM

Pathogenwatch is the public context provider. The pilot freezes metadata and
annotations once, downloads a bounded pool before SKA screening, and retains the same-ST,
divergence, recombination and temporal-signal checks. A failed dating gate is a
valid result: catalogue country composition and a topology/context report remain
useful without claiming a supported dated tree.

## Inputs and scope

Use a public Klebsiella pneumoniae ST147 catalogue (`organismId=573`, MLST field
`mlst`) containing *all* public QC-pass/fail records and verified detailed metadata.
Export cgLIN assignments for the same frozen UUID set with the
`scgMLST629_S` scheme. Keep raw search/details, export bytes, UTC retrieval time,
query/snapshot hashes and provider capability/version limitations.

Supply `focal.csv` with assemblies, valid dates, country/location, species, ST147,
origin `local`, and one explicit reference. Select three focal samples deliberately
from a documented related cgLIN prefix with strong source-ID/BioSample crosswalks;
record the selection rationale, exact IDs and assembly SHA-256s in the run README.
For a public demonstration, downloaded source genomes may act as a focal stand-in:
label them as such and exclude their accession aliases from context selection.
Do not describe them as a private survey. Real focal data require the user's
validated crosswalk/typing export; never copy a nearest public genome's cgLIN.

`focal_crosswalk.csv` contains `sample_id`, `source_genome_id` and, where available,
`biosample`. Check ambiguity/conflict/missing assignment output rather than
inferring cgLIN from ST. Root collection dates retain year/month/day or interval
precision; do not invent an exact day. Missing host/isolation source is retained.

## Laptop commands

Install the repository's Pixi environment (`pixi install`). From the checkout,
use the same Python environment for every stage. Keep the API key in protected
configuration or `PATHOGENWATCH_API_KEY`; it is never a command-line argument.

```bash
pixi run python scripts/st147_pilot.py --stage prepare \
  --focal pilot_inputs/focal.csv --catalogue pilot_inputs/st147_catalogue.json \
  --cglin-export pilot_inputs/st147_cglin.csv \
  --focal-crosswalk pilot_inputs/focal_crosswalk.csv \
  --output validation/st147_pathogenwatch/run --candidate-pool 8 --max-context 4 \
  --threads 2 --seed 20261009 --dry-run

# Bounded acquisition + SKA screening (up to 8 candidates, select up to 4).
pixi run python scripts/st147_pilot.py --stage prepare \
  --focal pilot_inputs/focal.csv --catalogue pilot_inputs/st147_catalogue.json \
  --cglin-export pilot_inputs/st147_cglin.csv \
  --focal-crosswalk pilot_inputs/focal_crosswalk.csv \
  --output validation/st147_pathogenwatch/run --candidate-pool 8 --max-context 4 \
  --threads 2 --seed 20261009

# Entirely offline selection + figure replay: no FASTA acquisition or native analyses.
pixi run python scripts/st147_pilot.py --stage replay \
  --context validation/st147_pathogenwatch/run/context \
  --output validation/st147_pathogenwatch/run

# Validate the analysis inputs/plan before opting into expensive native analysis.
pixi run python scripts/st147_pilot.py --stage analyse \
  --output validation/st147_pathogenwatch/run --threads 2 --dry-run

# Full native analysis, where supported, preserving default dating/divergence gates.
pixi run python scripts/st147_pilot.py --stage analyse \
  --output validation/st147_pathogenwatch/run --threads 2 --date-randomisations 100 \
  --public-focal-demonstration
```

Use `--public-focal-demonstration` only when the focal genomes are the public
stand-ins described here. It marks the report prominently as a demonstration;
omit it for a real focal cohort.

The minimum sample count is four, permitting this bounded demonstration rather
than changing temporal significance. The pilot uses p=0.05 and the existing
coherence/divergence gates. Record the resulting gate outcomes explicitly.
Native tools must be installed and supported on the execution platform; metadata,
selection and geography replay can be exercised separately on other laptops.

The driver records exact argument arrays, input SHA-256s, Python/platform, threads,
seed, return status and elapsed time in `prepare_provenance.json` or
`analyse_provenance.json`. Context provenance includes tool versions, candidate
and selected IDs, reasons, per-record download ledger, checksums, requests/bytes
where available, and rule-specific losses. The failure ledger must remain visible
when a partial download reduces the pool. Before release, verify resumability,
corrupt/interrupted downloads, exact source-ID-to-FASTA mapping and cache
invalidation with the offline test suites. A successful single download does not
prove the bulk endpoint is reliable.

## Same inputs on SLURM

Stage the complete pilot context, all referenced assemblies and checkout on the
cluster. Metadata/manifest assembly paths must resolve there; re-prepare from
frozen inputs on the cluster if paths differ. Install the same locked environment.
No heavy biological computation runs inside web requests.

```bash
export CHRONOCLADE_PILOT_ROOT=/absolute/path/to/staged/pilot
export CHRONOCLADE_PYTHON=/absolute/path/to/chronoclade/.pixi/envs/default/bin/python
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/st147_pilot.slurm
```

The script requests 4 CPUs, 16 GiB and 12 hours, launches the same analysis command
with 100 date permutations, and records the SLURM job log. Account/partition are
site-specific submission parameters. Preserve job ID, scheduler/resource settings,
software environment and input/output hashes. Cluster access is not required to
validate the script and laptop workflow; distinguish a documented runnable SLURM
path from an actually submitted/completed cluster run.

## Validation and offline reproduction

```bash
pixi run python -m pytest tests/test_pathogenwatch.py tests/test_cglin.py \
  tests/test_pathogenwatch_download.py tests/test_context_geography.py \
  tests/test_pathogenwatch_e2e.py -q
```

The end-to-end offline fixture freezes raw public responses, imports complete and
partial cgLIN, includes duplicate BioSamples, Unknown country, an undated record,
QC failure, singleton and a deliberately failed candidate download. It invokes
the real `prepare-context` CLI with download/SKA transports substituted, validates
manifest/combined metadata, checks full-public and selected geography, and replays
selected IDs/reasons plus every country table row without live credentials.
These fixtures are plumbing/contract tests, not evidence for biological
accuracy.

An opt-in upstream contract smoke test is separate from ordinary CI:

```bash
pixi run python scripts/st147_pilot.py --stage live-smoke --confirm-live \
  --output validation/st147_pathogenwatch/live_smoke
```

It checks public search and one UUID/details join without credentials. It does not
freeze all detailed metadata or verify download/cgLIN capability. Run protected
export/download smoke checks separately with a small known batch and record
expiry, polling, request count and partial failures. Normal CI never depends on
upstream uptime.

The offline replay checks the annotated snapshot hash, reconstructs selection from
stored SKA distances and settings, and compares exact selected source IDs/reasons
and source table counts. Changing a snapshot, cgLIN scheme/version, selection
settings or focal assembly requires affected preparation/selection stages to be
rerun, with a changed selection fingerprint. Never silently reuse another run's
manifest. Regeneration of country figures does not require phylogeny.

Before release, hand-check a subset of raw-export groups/countries/focal memberships;
reconcile eligible genome/sample-unit counts with duplicate audit, valid groups and
unresolved coverage. Inspect `context/context_geography/index.html`, the integrated
report, exported SVG/PNG and supporting archive at readable size. The full public
catalogue denominator includes undated QC-pass records and is never replaced by
the selected tree tips. Keep separate public, selected and focal denominators.

The figures measure composition of available sequenced records, not prevalence,
incidence, migration direction or transmission. Submission/surveillance bias,
incomplete or provisional cgLIN, organism-specific capability, metadata conflicts
and date uncertainty remain material limitations. Country is not country of
acquisition. A same-ST catalogue is not the full global lineage across all STs.

## Committed public freeze for offline figures

The permitted compact public fixture retains all 5,647 QC-passing deduplicated
units, their source UUIDs, sample-unit/raw-record counts, normalised country/date/QC
fields, full cgLIN scheme/version/status/prefix annotations, three verified public
focal stand-ins and retrieval/source hashes. Raw detailed responses and assemblies
are excluded; this file is a geography/coverage freeze, not a replacement for the
complete metadata source envelope.

```bash
pixi run python scripts/st147_pilot.py --stage figures \
  --frozen-fixture validation/st147_pathogenwatch/frozen_catalogue.json.gz \
  --output validation/st147_pathogenwatch/offline_figures
```

This validates the compact payload hash and creates figures/tables/HTML with no
network, assemblies or phylogenetic analyses. The corresponding real-freeze test
hand-checks the focal group's country subset and denominators, all-depth assignment
coverage, raw/sample-unit reconciliation, zero unresolved sample identity and three
focal memberships/overlaps. The source's immutable raw search/details and annotated
catalogue hashes remain embedded for audit.

Metadata and manifest assembly/catalogue paths are written relative to their
files. For SLURM analysis, stage the complete pilot directory, including focal
assemblies, context assemblies and annotated catalogue, preserving that layout.
The same frozen directory can move between laptop and cluster without rewriting
absolute laptop paths; a regression fixture verifies the moved inputs.
