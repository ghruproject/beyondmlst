# Experimental ESM2 protein embeddings

The optional ESM2 module embeds each unique amino-acid sequence in a supplied
protein FASTA once and saves vectors with a complete source-record mapping.
An optional prepared-dataset and sample–locus mapping route adds a descriptive
date-distance report alongside a conventional cgMLST NJ baseline. It can reuse
saved embeddings without loading the model. The prepared-data genome route resolves a supplied scope-checked DNA catalogue,
translates CDS, computes locus-preserving genome distances and exports the same
validated selection/ensemble contracts as cgMLST. Live catalogue acquisition and
corrected-SNP biological validation remain separate work in the
[tool suite design](tool-suite-design.md).

The output is experimental. Protein embedding distances have different units
from allele or SNP differences. Runtime benchmarks do not demonstrate neighbour
or cluster preservation, and neither model is a validated scientific default.
The command starts with 8M as a small operational setting; model selection is explicit.

## Installation and local execution

The base ChronoClade install and CLI help do not import PyTorch or fair-esm, acquire
model weights, or require a GPU. Install the optional dependencies in an isolated
Python 3.11 environment to avoid changing the all-tools Pixi environment:

```bash
uv venv --python 3.11 .venv-esm2
uv pip install --python .venv-esm2/bin/python --editable '.[esm2]'
.venv-esm2/bin/chronoclade esm2 --help
```

Supported upstream model names are `esm2_t6_8M_UR50D` (8M, 320 dimensions, layer 6)
and `esm2_t12_35M_UR50D` (35M, 480 dimensions, layer 12).
The aliases `8M` and `35M` select these exact architectures. Extraction uses
float32 and the mean of the final layer's residue representations. BOS, EOS and
padding tokens are excluded. Contact prediction and ESMFold are not used.

Supply explicitly acquired local checkpoint weights for offline execution:

```bash
.venv-esm2/bin/chronoclade esm2 proteins.faa \
  --out esm2-8m --model 8M --device cpu \
  --checkpoint /path/to/esm2_t6_8M_UR50D.pt \
  --checkpoint-sha256 YOUR_64_CHARACTER_CHECKPOINT_SHA256 \
  --token-budget 4096 --cache-dir embedding-cache
```

The loader reads the named file directly and never invokes the hub loader on
this path. It validates the checkpoint architecture and any supplied SHA256.
It does not require the contact-regression sidecar, because contact predictions
are not requested. Official checkpoints contain Python configuration objects;
open only weights obtained from a trusted source. A file hash pins the acquired
file but is not independently a statement about its authenticity.

An explicit `--allow-download` permits retrieval of the selected upstream model
into PyTorch's checkpoint directory when no local checkpoint is supplied:

```bash
.venv-esm2/bin/chronoclade esm2 proteins.faa \
  --out esm2-35m --model 35M --device auto --allow-download
```

Weights are never downloaded during base installation, CLI startup or an offline
checkpoint run. The download branch uses the official model URL and is opt-in.
The manifest records the acquired file's SHA256 and the runtime versions.

## Inputs, batching and caching

Input is an uncompressed UTF-8 protein FASTA with unique record IDs. Full headers
are preserved, so caller-provided genome and locus identifiers remain inspectable;
they are not automatically parsed into an invented genome/locus contract.
Lower-case residues are normalised to upper case. Standard residues plus ESM's
`X`, `B`, `Z`, `U` and `O` are accepted. Empty sequences, stops, alignment gaps,
internal spaces, duplicate record IDs and malformed FASTA fail clearly. Input
validation is a protein alphabet check, not a DNA translation or coding-sequence
validation service.

The default and maximum supported protein length is 1022 residues. Longer
sequences fail; extraction never silently truncates them. A lower `--max-length`
can restrict a run. The library's `invalid_policy="exclude"` explicitly retains
invalid-sequence exclusions in the manifest while embedding the remaining valid
records. Structural FASTA errors and duplicate IDs always fail.

Batches are ordered by length and bounded by padded token count, including BOS
and EOS: `batch size × (longest protein length + 2) <= token budget`. A single
protein exceeding the budget fails before model loading. This bounds tokens,
not transformer memory: long sequences still incur quadratic attention cost.

Identical normalised proteins share one SHA256 identifier and one vector.
Every retained source record maps back to that identifier, including different
DNA alleles that translate to the same protein when supplied by an upstream caller.
No synonymous DNA distinction can be recovered from an amino-acid embedding.

Optional cache entries include the protein ID, vector checksum and full
extraction provenance: model/checkpoint hash, representation layer, pooling,
precision, runtime versions, actual device and token budget. Changed provenance
uses a separate cache namespace. Reuse validates identity, dimensions, dtype,
finiteness and checksums without enabling pickle loading. Corrupt entries are
recomputed and counted in the audit. Complete extraction is required before new
cache entries are published; CPU fallback cannot publish partially computed GPU
vectors under CPU provenance.

Devices are `cpu`, `mps`, `cuda` and `auto`. Auto chooses CUDA, then available
Apple MPS, then CPU. If accelerator allocation or extraction fails, auto retries
on CPU and records the reason and actual device. Explicit unavailable devices
or failed explicit inference fail clearly. Actual MPS compatibility must be
demonstrated by a completed run, rather than inferred from an availability check.

## Saved artifacts and Python API

Each successful run writes:

- `embeddings.<content-hash>.npz`: `protein_ids`, float32 `embeddings` and residue
  `lengths`. Row order matches `protein_ids` in the manifest.
- `embeddings.json`: schema/software versions, input checksum and scope, exact
  IDs, all source mappings, exclusions, model/runtime/device provenance, batch
  and cache counts, timings, and the vector artifact's relative path and checksum.

Artifacts use atomic file replacement. The numerical artifact has an immutable
name; the completion manifest is replaced last. An interrupted rerun therefore
does not invalidate the preceding completed manifest. Previous numerical
artifacts may remain in the output directory after a changed rerun.

```python
from pathlib import Path
from chronoclade.esm2 import run_embeddings

result = run_embeddings(
    Path("proteins.faa"), Path("esm2-results"),
    model="8M", checkpoint=Path("esm2_t6_8M_UR50D.pt"),
    device="cpu", token_budget=4096, cache_dir=Path("embedding-cache"),
)
print(result.manifest_path)
print(result.vectors_path)
```

## Prepared-dataset genome exploration

See [the genome workflow](esm2-genomes.md) for independent prepared-data analysis,
real DNA translation, fixed-panel distances, shared metadata ordination controls,
nearest neighbours, experimental groups and selection ensembles.

## Prepared-dataset date-distance report

Supply the prepared `dataset.json` and a CSV with exactly named
`sample_id,locus,record_id` columns. Each row links a stable prepared sample ID and
a locus in its compatible cgMLST catalogue to a source record in the embedding
manifest. Record IDs identify the supplied FASTA proteins; they are not guessed
from headers. Mapping identities, scheme compatibility, vector checksums and
model provenance are validated before reporting. Eligible mapped samples must belong
to one species; mixed-species comparisons fail rather than crossing species with
a single embedding reference.

To create embeddings and the report together:

```bash
.venv-esm2/bin/chronoclade esm2 proteins.faa \
  --dataset prepared/dataset.json --sample-loci sample-loci.csv \
  --out esm2-temporal --model 8M --device cpu \
  --checkpoint /path/to/esm2_t6_8M_UR50D.pt
```

To analyse existing vectors without importing PyTorch or fair-esm:

```bash
chronoclade esm2 --embeddings-manifest esm2-8m/embeddings.json \
  --dataset prepared/dataset.json --sample-loci sample-loci.csv \
  --out esm2-temporal --reference-sample YOUR_STABLE_SAMPLE_ID \
  --panel-locus locus_A --panel-locus locus_B
```

The reference is the requested real sample or the lexicographically first
eligible stable sample ID. It is chosen independently of collection dates. The
protein panel is explicit with repeated `--panel-locus` options, or uses the
mapped-locus intersection. A sample must have a valid mapped protein at every
panel locus to receive a distance; missing embeddings are never zero vectors.
The numerical artifact records panel IDs, missing loci, exclusions and per-locus
distances. An incomplete or empty panel fails clearly.

`index.html` uses the existing offline ChronoClade full-report stylesheet. It has
two named measurements with shared location colours and input/context markers:

- **cgMLST rooted NJ distance versus collection date:** conventional categorical
  mismatch fractions on one fixed set of loci callable in every profile in each
  compatible cohort. The NJ tree is rooted at the same real reference sample;
  incompatible cohorts remain separate. Negative NJ branches are retained and
  counted. Slope units are **cgMLST mismatch fraction/year**.
- **Protein distance to reference versus collection date:** the arithmetic mean
  of per-locus cosine distances to the fixed reference on the frozen protein
  panel. Slope units are **mean locus cosine distance/year**.

Both views use unweighted ordinary least squares with collection-date interval
midpoints, retaining day, month, year and interval bounds as horizontal bars.
Bars describe date precision, not statistical confidence intervals. Each figure
includes residuals in its own distance units. Finite distances without eligible
dates remain in an undated panel; all sample rows and exclusions stay in the CSV
and JSON. Fewer than three distinct date midpoints leaves the fit unavailable.
Constant distances have slope zero and undefined correlation/R², shown explicitly.

The report records prepared, mapped, panel-eligible and dated sample counts, the
reference/root, callable loci, scheme/database versions and actual model and
runtime provenance. Downloads include SVG/PNG figures, per-panel CSV files,
`temporal_points.csv` and complete `temporal_diagnostics.json` evidence.
`temporal_analysis.json` fingerprints the inputs, parameters and completed report.

These diagnostics are exploratory, with no formal clock test, confidence
intervals or p-values. Shared ancestry and sampling structure can produce date
association. Protein embeddings cannot distinguish synonymous DNA changes that
produce identical proteins. Neither the protein distance nor its fitted slope is
a substitution rate, and no ancestry or TMRCA is inferred from embeddings.
Supported dating remains the separate assembly-based tree/time route.

## Separate runtime benchmark

The benchmark module compares the exact same frozen protein FASTA across models
and devices. It warms up inference, then records fresh cache-disabled repetitions;
cache population and validated reuse are measured separately.

```bash
.venv-esm2/bin/python -m chronoclade.esm2.benchmark proteins.faa \
  --out esm2-runtime-comparison \
  --checkpoint-8m /weights/esm2_t6_8M_UR50D.pt \
  --checkpoint-35m /weights/esm2_t12_35M_UR50D.pt \
  --devices cpu mps --token-budget 4096 --repeats 3
```

Choose only devices available on the benchmark host. Unavailable or failed
explicit combinations are recorded as failures and result in a partial report.
`benchmark.json` retains raw fresh-inference timings, median residue throughput,
checkpoint/model-load time, runtime and weight provenance, length distribution,
batch/token counts, process peak RSS and accelerator memory where measurable.
Accelerators are synchronised before timing boundaries. Process peak RSS is
cumulative within the process; MPS allocation snapshots are not peak-memory
measurements. Imported libraries and operating-system file caches remain warm
for later model loads, so these load times are not a controlled cold-disk study.

Use real protein alleles, close variants and repeated proteins derived from
synonymous DNA alleles for later validation. Runtime results alone cannot establish
the scientific adequacy of 8M or 35M. Biological neighbour/cluster preservation,
missing-locus robustness and genome selection coverage require separate studies.

### Measured laptop starting setting

On 10 October 2026, both models ran successfully on an Apple M1 Pro with 32 GiB
memory using CPU and native MPS. The frozen fixture contains 32 unique proteins
translated from public Pasteur `scgMLST629_S` allele-1 records, 102–748 residues
long. Bacterial genetic code 11 and CDS initiation rules were used. It is a
runtime fixture, not a panel of 32 genomes or a clustering validation dataset.

Median warm extraction times from three repetitions, after a full-fixture
warmup, with float32, four CPU threads and a 4096 padded-token budget:

| Model | CPU | Apple MPS GPU |
| --- | ---: | ---: |
| 8M | 1.395 s | 0.501 s |
| 35M | 3.370 s | 1.137 s |

These extraction timings exclude checkpoint download and model loading. The
separate model loads, observed process memory, sequence provenance and numerical
checks are retained in the
[frozen benchmark evidence](https://github.com/ghruproject/chronoclade/tree/main/validation/esm2).
Both models produced finite vectors on both devices, and cache reuse retained
identical vectors. This supports **8M as the initial operational default** for
speed and memory; 35M remains explicit. The relative clustering usefulness of
the models remains unmeasured, and timings should not be extrapolated directly
to a 10,000-genome analysis without counting unique proteins and their lengths.

## Offline SLURM example

Stage the environment, protein FASTA and trusted checkpoints on a login or
staging node before submission. Record the checkpoint hash and use it in the job.
Adapt partition, account, paths, memory and wall time to the site's configuration.
The following CPU example uses no GPU allocation:

```bash
#!/bin/bash
#SBATCH --job-name=chronoclade-esm2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Add site-specific --partition and --account directives here.
set -euo pipefail

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

RUNTIME=/project/chronoclade/.venv-esm2/bin/chronoclade
PROTEINS=/project/input/proteins.faa
WEIGHTS=/project/weights/esm2_t6_8M_UR50D.pt
CHECKPOINT_SHA256=YOUR_64_CHARACTER_CHECKPOINT_SHA256
JOB_OUTPUT=/project/results/esm2-${SLURM_JOB_ID}
JOB_CACHE=/scratch/${USER}/esm2-cache-${SLURM_JOB_ID}

"$RUNTIME" esm2 "$PROTEINS" --out "$JOB_OUTPUT" \
  --model 8M --device cpu --checkpoint "$WEIGHTS" \
  --checkpoint-sha256 "$CHECKPOINT_SHA256" \
  --token-budget 4096 --cache-dir "$JOB_CACHE"
```

For CUDA, provision a compatible GPU runtime, request one GPU using the site's
resource directive, and select `--device cuda`. An array should shard unique
protein sequences through explicit input manifests, give workers separate output
and cache namespaces, and validate completeness and hashes before assembling
shards. This package does not yet provide automatic shard assembly. This script
is documentation; it does not claim that a real cluster run has been tested.

Upstream implementation references:
[ESM extraction examples](https://github.com/facebookresearch/esm#extracting-embeddings),
[checkpoint loading](https://github.com/facebookresearch/esm/blob/main/esm/pretrained.py)
and [batch conversion](https://github.com/facebookresearch/esm/blob/main/esm/data.py).
