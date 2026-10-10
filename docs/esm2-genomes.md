# Prepared-data ESM2 genome exploration

The ESM2 route is an independent alternative to `cgmlst` after `prepare`. It uses
actual allele coding sequences and keeps locus identity throughout genome
comparison. Its selections have the same validated contract consumed by `tree`.
The route remains experimental: numerical correctness and model execution are
exercised, but corrected-SNP within-lineage biological adequacy is not established.

```bash
chronoclade esm2 prepared/dataset.json --out esm2-results \
  --allele-catalogue frozen-dna/catalogue.json \
  --model 8M --device mps --checkpoint /weights/esm2_t6_8M_UR50D.pt \
  --selection-size 50 --selection-runs 5 --nearest-per-input 3
```

The catalogue contract is documented in [frozen DNA inputs](esm2-sequence-input.md).
A compatible known database fingerprint is required. Catalogue DNA is translated
with bacterial genetic code 11, CDS initiation semantics and terminal-stop
validation. Missing sequences and invalid CDS are audited; alleles and hashes are
lookup keys and cannot reconstruct DNA. `--terminal-stop allow_absent` explicitly
permits stop-less caller outputs; no gaps, frames or internal stops are repaired.
Weights are local unless `--allow-download` is explicitly supplied. CPU, MPS and
CUDA are supported; explicit unavailable devices fail and automatic fallback
records the actual device. Model libraries remain lazy and are unnecessary when
analysing already saved vectors.

Saved-vector execution, without importing PyTorch or fair-esm:

```bash
chronoclade esm2 prepared/dataset.json --out esm2-reused \
  --embeddings-manifest prior/embeddings.json --sample-loci sample_loci.csv \
  --panel-locus locus_A --panel-locus locus_B
```

The original protein FASTA extraction command and optional `--dataset` temporal
report remain available. Use `--genome-analysis --dataset prepared/dataset.json`
when supplying the dataset as an option instead of positionally.

## Distances, missing loci and groups

The route uses the same shared, explicit disjoint lineage partition policy as
cgMLST: `--lin-level 5`, `6` or `7`, or a scheme-specific `--hiercc-level HC10`.
Every matching context record already present in the frozen bundle is considered;
the command does not claim fresh or globally complete public retrieval. Unresolved
inputs and unmatched contexts remain in `partitions.json` and the index.
Each resolved block requires one compatible profile matrix and catalogue scope.

Unique proteins are embedded once per block; a cache permits reuse across blocks.
Genome distance is the equal-weight mean of locus-specific cosine distances on
one fixed panel. The panel uses the mapped-locus intersection across samples with
any mappings, or explicit repeated `--panel-locus` values. Every scheme locus
remains in the saved protein-ID table and mapped mask. Incomplete explicit-panel
samples have no distance; missing vectors are never zero vectors. Panel omission,
variable protein loci and all exclusions remain inspectable.

Ordination uses `sqrt(2 × mean locus cosine distance)`: the Euclidean chord
geometry of locus-preserving unit protein vectors. Axis labels identify protein
chord units. Country, date, role and group views use the same shared offline
ordination viewer as cgMLST. All eligible genomes retain their identities.
`--neighbour-count` retains every tie at its requested boundary. Groups are
connected components at `--group-threshold` in embedding-distance units. The
default zero threshold groups identical protein panels; it is not a calibrated
clonal-group threshold or a substitute for conventional phylogeny.

The current experimental dense genome backend is explicitly limited to 1,500
eligible genomes. It fails at larger sizes rather than silently subsampling.
This is separate from improvements to the conventional large-profile backend.

## Selections and reports

The input genomes, pinned exact prepared IDs (`--pin`) and required nearest
context ties are mandatory. Shared selection policies retain budget overruns,
per-ID reasons, metadata/group strata, fixed quotas and reproducible alternative
seeds. Selection distance evidence includes checkpoint, vector and fixed-panel
fingerprints, with ESM2 units explicitly recorded. A missing input or pinned
protein panel produces an unavailable selection with exact IDs in the report;
mandatory genomes are never silently omitted. Ineligible context genomes are
excluded from the selection pool with an explicit mask and coverage audit.

Each block saves:

- `genome_analysis.json`: exact IDs, panel, model/checkpoint/runtime evidence,
  nearest neighbours, groups, exclusions and biological comparison status.
- `genome_vectors.npz`: pair distances, chord ordination, all sample/locus protein
  IDs and the complete mapped-locus mask.
- `pairwise_distances.csv`: ESM2 and allele comparisons on the same fixed panel.
  Catalogue-driven runs also report mean per-locus DNA mismatch fractions for
  equal-length CDS. Length-changing pairs are marked unavailable; no unrecorded
  alignment is invented.
- `selection_bundle/ensemble.json` and its checksummed individual selections,
  validated by the same loader as cgMLST selections and accepted by `tree`.
- `index.html`, shared metadata ordination views, and labelled baseline reports.

By default, the baseline report computes the conventional cgMLST NJ guide tree
and approved location-state network from categorical alleles. The collection-date
report separately compares cgMLST root-to-tip and protein-to-reference distances,
with interval bars and residuals. Both reuse existing report components and style.
The baseline tree and network have their own membership and missingness evidence;
embedding groups and coordinates never drive biological reconstruction.
`--no-baseline` omits those conventional/date reports.

The stage publishes `esm2.json` last, with artifact checksums, source hashes,
parameters, block statuses and report links. It uses a fresh immutable output
namespace; interruption leaves previous outputs intact. Selection bundles retain
portable relative references when moved together with the prepared dataset.

## Executed real-sequence validation

On 10 October 2026, 8M and 35M both completed actual CPU and Apple MPS inference
on two frozen public Pasteur bacterial CDS and explicit single-codon synonymous
and nonsynonymous controls. Six simulated genome combinations included one missing
context locus; five genomes were eligible on a fixed two-locus panel. The fixtures
use controlled metadata and lineage assignments and are not a sampled population.

Across all four combinations, the synonymous control retained two allele
mismatches and nonzero DNA mismatch distance while producing exactly zero ESM2
distance. The missing context remained excluded and audited; input and pinned
synonymous IDs were retained in validated selection ensembles. The full 8M CPU
report included the conventional baseline and collection-date panels. Country and
date ordination controls were exercised in the rendered browser report.

The complete [prepared fixture](esm2-genome-fixture/prepared/dataset.json) and
[DNA catalogue](esm2-genome-fixture/catalogue/catalogue.json) are frozen with these
controls, so the genome workflow can be repeated using local checkpoint weights.
Run the command above with those paths and repeated `--panel-locus hfq_S` and
`--panel-locus KP1_2532_S` to retain the incomplete-context test.

[Saved validation evidence](esm2-genome-validation.json) contains checkpoint and
runtime provenance, timings and allele neighbour agreement. These small-fixture
inference times include first-use accelerator overhead and are not a representative
throughput comparison; use the larger runtime benchmark documented in
[ESM2 embeddings](esm2.md). Agreement on these controls does not validate either
model for within-lineage selection. Corrected SNP neighbour/cluster validation and
missing-locus sensitivity on representative natural genomes remain required before
adopting a scientific default. Embeddings cannot infer synonymous DNA change,
ancestral states, mutation rates or dates of common ancestors.
