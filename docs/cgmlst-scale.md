# Full-profile large-block cgMLST analysis

Independent `cgmlst` blocks containing more than 1,500 records use the scale
backend. Every supplied identifier stays in the report and binary evidence.
Every eligible profile enters a complete-comparability cohort, its full NJ tree,
ordination and metadata summaries. No sample-size cap, representative-tree
substitution or silent sampling is used to claim complete analysis.

## Installation and use

```bash
pixi install -e cgmlst
pixi run -e cgmlst chronoclade setup-scale
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --output cgmlst
```

The isolated `cgmlst` environment includes SciPy. On Linux x86_64 it includes the
Bioconda RapidNJ package. On supported macOS hosts `setup-scale` builds the
[pinned upstream source](https://github.com/somme89/rapidNJ/tree/ed2d36e219d9db16778b941b5054c0fd021b528a)
without modifying its algorithms and checks a real three-tip distance tree.
The source archive SHA-256 is
`57ea662a3589459528beaa96de011bd1b389630484d9e3215a942bdb19518f21`.
A build log and installation receipt record the source revision, platform,
compiler command and executable checksum. The default install is
`~/.local/share/chronoclade/scale/bin/rapidnj`.

The upstream source uses x86 SSE. Apple Silicon builds an x86_64 executable and
requires an already installed Rosetta runtime. Linux ARM source installation
fails explicitly; use a separately validated executable. Select a backend with
`CHRONOCLADE_RAPIDNJ`, put `rapidnj` on `PATH`, or use the verified local install.
Missing or failed RapidNJ never triggers a replacement tree or sampled analysis.

The CLI's automatic bootstrap default is 30 for small blocks and zero for large
blocks, recorded per block. The scale backend currently does not implement locus
bootstrap grouping: an explicit positive request fails before writing the
analysis, with an actionable explanation. It does not omit requested replicates.

## Distances and comparability

Each locus is encoded independently as categorical labels. Batched compiled
SciPy Hamming kernels count all encoded mismatches and the presence-mask XOR.
Subtracting the XOR excludes exactly the loci called in only one profile;
call totals and XOR give the jointly called denominator. Distances are mismatch
counts divided by jointly called counts. Numeric allele identifiers are never
ordinary quantitative features.

Species, scheme/version, database fingerprints and locus universes retain the
existing compatibility rules. Observed export locus unions remain explicitly
labelled when canonical completeness is unknown. Missing calls are not matches.
Insufficient-overlap and incompatible comparisons are NaN, never zero-imputed.
Self-ineligible profiles remain in exclusions and metadata. Deterministic greedy
cohorts contain only complete pairwise comparability; valid cross-cohort pair
evidence remains available in the global matrices.

The disk-backed evidence directory contains four square NumPy matrices:

| File | Type | Meaning |
| --- | --- | --- |
| `distance.npy` | float64 | mismatch fraction; unavailable pairs are NaN |
| `allele_differences.npy` | uint32 | observed called mismatch count |
| `shared_called_loci.npy` | uint32 | jointly called denominator |
| `scheme_loci.npy` | uint32 | declared or explicitly observed locus-universe size |

`distance_evidence.json` records sample order, scope, cohort membership, unavailable
pair counts, paths and checksums. The loader verifies identifiers, scope, file
checksums, matrix shapes/types, symmetry, finite ranges, overlap, comparable-pair
counts and agreement between normalized distances and callable counts. Selections
use lazy matrix pair lookup rather than quadratic Python dictionaries; their
saved evidence bundle includes the binary files and a portable manifest.

All nearest public comparisons rank by **raw called allele differences** among
pairs meeting the overlap rule. Raw ties are retained along with their actual
normalized distance and denominator. Such denominators must be inspected when
comparing incomplete profiles.

Large analyses do not automatically write every pair to CSV. To export explicitly
requested rows:

```python
from chronoclade.matrix_distances import BinaryDistanceEvidence

evidence = BinaryDistanceEvidence("binary_distances/distance_evidence.json")
evidence.export_pairs("requested.csv", pairs=[("input-ID", "public-ID")])
# Or all pairs among an explicitly supplied list of IDs:
evidence.export_pairs("selected.csv", sample_ids=["a", "b", "c"])
```

Unavailable requested comparisons remain marked unavailable with their preserved
callable-count evidence where available. They do not acquire an invented distance.

## Tree, ordination, groups and country network

Each cohort also stores a square distance matrix and condensed vector. Complete
linkage grouping uses compiled SciPy routines and reports threshold sensitivity.
Input ordering fixes its tie behaviour, although an exact tied merge may differ
from the legacy Python implementation. The groups remain descriptive, without
claiming outbreak, transmission or biological calibration.

RapidNJ consumes streamed PHYLIP distances with safe temporary tip names, restores
original IDs and verifies the complete output tip set. The Newick is downloadable.
Input values are written with 17 significant digits; upstream RapidNJ uses float32
internally. Negative branches are retained and flagged. The temporary text matrix
is removed after a successful tree, while the binary evidence remains.

PCoA uses ARPACK leading eigenpairs of the exact centered squared-distance operator.
Each matrix operation runs in bounded row batches; there is no dense full Gram
matrix or full eigendecomposition. All cohort genomes receive coordinates. The
saved solver tolerance is explicit. Total positive inertia and the complete
negative spectrum are **not computed**, and no fabricated variance percentage is
shown. Iterative convergence failure remains a failure.

The full tree is rooted at the lexicographically first input tip, or first eligible
tip in a cohort without inputs. Root provenance and sensitivity are recorded.
Exploratory date diagnostics use declared allele-distance units on a fixed fully
called panel. When callable denominators vary, the mismatch fraction is multiplied
by the recorded locus-universe size and labelled **scheme-equivalent** allele
units per year; it is not an observed mutation or allele-change rate. Date
intervals, the fitted line, unsupported slopes and exclusions remain visible.

Country reconstruction uses the existing full-tree state-change network, exact
optimal-history pair ranges, root sensitivity, centrality metrics and approved
styling. Indexed exact tree distances accelerate root selection without changing
its answers. The standalone country-coloured tree lazily expands clades, retaining
all IDs and repeated profiles. It uses the same representative ancestral country
states and palette as the network. No enormous static all-tip SVG is generated.

## Measured synthetic resource evidence

Runs on 10 October 2026 used macOS 14.6.1/Apple Silicon, the verified x86_64 RapidNJ
source build under Rosetta, SciPy 1.17.1, 128-profile batches and a 2,048 MB RapidNJ
memory allowance. Each labelled synthetic workload used 629 loci, 20 categorical
templates, three altered loci per profile, ten inputs and eight collection years.
These are computational checks, not biological validation or SNP-pipeline timings.

| Profiles | Known countries | Full analysis | Distances | Full network and uncertainty | Peak Python RSS | Peak RapidNJ RSS | Output disk |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1,500 | 4 + unknown | 8.16 s | 3.20 s | 1.21 s | 385 MB | 25.7 MB | 80.9 MB |
| 1,500 | 20 + unknown, independently assigned | 12.47 s | 3.22 s | 5.36 s | 388 MB | 25.6 MB | 83.0 MB |
| 10,000 | 4 + unknown | 99.74 s | 36.41 s | 5.55 s | 3.45 GB | 862 MB | 3.25 GB |

Every run verified that available-profile counts, full NJ tip counts and PCoA row
counts equal the requested size. The 10,000-profile run retained all 49,995,000
pair comparisons. Network measurements include exact fixed-tree uncertainty and
five candidate roots. The 10,000-profile workload's modest country cardinality
and structured distances do not demonstrate arbitrary high-country or random
singleton-group performance. Broader biological and cluster benchmarks remain
necessary before stronger resource or scientific claims.

Evidence: [1,500 profiles](cgmlst-scale-benchmark-1500.json),
[1,500 profiles / 20 countries](cgmlst-scale-benchmark-1500-20countries.json),
[10,000 profiles](cgmlst-scale-benchmark-10000.json).

Reproduce the labelled workloads:

```bash
pixi run python scripts/benchmark_cgmlst_scale.py --samples 1500 --output /tmp/scale-1500
pixi run python scripts/benchmark_cgmlst_scale.py --samples 10000 --output /tmp/scale-10000
pixi run python scripts/benchmark_cgmlst_scale.py --samples 1500 --countries 21 --country-independent --output /tmp/scale-1500-20countries
```

Exact pair computation and storage remain quadratic. Global evidence needs
approximately `20 × n²` bytes; one full cohort adds `12 × n²` bytes for square and
condensed distance storage. A portable selection copy adds another `20 × n²`
bytes. The main 10,000-profile analysis therefore requires about 3.2 GB of matrix
storage before selection copies, temporary PHYLIP, reports and filesystem overhead.
PCoA batches are bounded, while compiled complete linkage and RapidNJ also require
quadratic working memory. The observed process peaks do not constitute a hard
whole-pipeline memory cap. Tens of thousands may require substantial cluster RAM
and scratch storage; no larger-size completion or scheduler allocation is claimed
from these runs.
