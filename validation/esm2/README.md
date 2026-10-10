# ESM2 laptop runtime evidence

Measured on 10 October 2026 on an Apple M1 Pro with 32 GiB memory. This tests
protein embedding inference and cache correctness, not genome clustering quality.

`proteins.fasta` contains 32 unique proteins translated from public Pasteur
scgMLST629_S allele-1 records using bacterial genetic code 11 and CDS initiation
rules. `fixture-manifest.json` records the source URLs, DNA hashes and translation
provenance. `correctness.fasta` additionally contains an exact duplicate and a
labelled synthetic substitution; it is not a panel of clinical isolates.

`benchmark-summary.json` contains all four model/device runs, raw repeated timings,
versions, checkpoint hashes, observed memory, CPU/MPS numerical comparisons and
actual package API/cache checks. Local paths in the captured evidence identify
the original run directory. Model weights, vector caches and virtual environments
are not committed.

| Model | CPU median | MPS median |
| --- | ---: | ---: |
| 8M | 1.395 s | 0.501 s |
| 35M | 3.370 s | 1.137 s |

These are warm extraction times for 32 proteins totalling 8,910 residues, after
one complete warmup and across three timed repetitions. Both models use float32,
residue-only mean pooling, four CPU threads and a 4096 padded-token budget.
Checkpoint acquisition and model loading are reported separately. The actual
embedding engine ran successfully on CPU and native MPS, with finite vectors,
correct deduplication and bit-identical vector-cache reuse.

8M is the initial operational default on this host. Neither model has yet been
validated for within-LIN-group genome clustering or representative selection.
CUDA, SLURM and a remote cluster were not tested.

For reproduction, see `scripts/benchmark_esm2.py` and the optional module's
`python -m chronoclade.esm2.benchmark` entry point. Install the optional ESM2
environment first; the standalone script also needs `psutil`. Use the recorded
local checkpoints and this frozen FASTA to avoid network activity during timing.
