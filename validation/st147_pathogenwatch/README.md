# ST147 Pathogenwatch release pilot

The executable laptop and SLURM recipe is in [the pilot guide](../../docs/pathogenwatch-pilot.md).
Use `scripts/st147_pilot.py` for bounded preparation, analysis and completely offline
selection/country-table replay. Ordinary CI generates compact synthetic public
fixtures in `tests/test_pathogenwatch_e2e.py`; no private material or credentials are
committed.

The public metadata catalogue frozen on **9 October 2026** contains 7,807 genome
records, 5,721 accession-deduplicated units before QC, and 7,711 explicitly QC-passing
records. It measures 7,734 records with country/73 Unknown, and date precision of
4,596 day, 46 month, 2,871 year, 260 missing and 34 invalid. These are raw record counts,
not patients or infections. The complete source envelope content hash after canonical Myanmar/UAE alias
normalisation and complete SAMN/SAMEA BioSample parsing is
`bcffd21c69089623103170b0e770f6852db3178be0a07a2d029176f94eb4b8d7`.
Raw responses are unchanged; re-normalisation removed all 39 apparent country
conflicts because those values differed only by aliases. Every raw record and every
QC-passing deduplicated unit has a recognised BioSample after SAMEA parsing. Use the saved run
manifest/provenance to verify the actual input file hashes.

The full annotated **QC-passing** public catalogue deduplicates to 5,647 sample units.
Measured coverage and group sizes are:

| Full-prefix depth | Assigned units | Unresolved units | Groups | Singleton groups | Largest group |
|---|---:|---:|---:|---:|---:|
| 5 | 5,558 | 89 | 104 | 76 | 4,371 |
| 6 | 5,546 | 101 | 240 | 169 | 2,558 |
| 7 | 5,546 | 101 | 635 | 441 | 1,767 |

Duplicate biological-sample assignment conflicts are conservatively unresolved at
every depth: 88 conflicted units plus 1 partial unit at depth 5, and 13 partial
units at depths 6/7. Raw genome counts reconcile as 7,711 records across 5,647
units.

Depth 5 is the initial overview because it has the highest assignment coverage and
fewest rows. Deeper views remain available in expandable sections and source tables.
Median group size is one at each depth: singleton-rich available sampling should not
be hidden. These prefix depths are not pairwise SNP cutoffs or transmission clusters.
Scheme `scgMLST629_S`; upstream nomenclature version is unavailable and recorded
`unknown`. Figures must be regenerated after any catalogue/annotation change.

A nominated public focal stand-in set has a common resolved depth-7 prefix
`[0,0,197,0,4,1,0]` and strong BioSample evidence:

| Source genome ID | BioSample | Collection date | Country |
|---|---|---|---|
| aLZnApmwKX44bLixaxccQn | SAMD00193192 | 2016-01-01 | Myanmar (Burma) |
| 2iNpme3TkKQVv7MFwugPF9 | SAMEA9453699 | 2019 | Italy |
| 8HFG63i6tfFmeqns6auhFF | SAMEA122094000 | 2022 | India |

They were selected from QC-passing dated public entries with valid shared depth-7
annotation, different years and strong accession links. Each has exactly one public
source genome record and was checked through the actual focal resolver against all
7,807 annotated records: complete assignments at prefix depths 5, 6 and 7. An earlier
raw-source nomination was rejected when the focal resolver exposed contradictory
assignments among BioSample aliases; no contradictory code was silently borrowed. This is a public demonstration,
not a private focal survey. Their public records remain once in the full-public
composition; their aliases must be excluded from context selection. The worked run
must verify downloaded assemblies and retain source/checksum crosswalks before
calling the focal acquisition completed.

Preserve the full raw snapshot and cgLIN export outside Git, along with file hashes,
exact queries, retrieval UTC, capability/version limitations, download ledger and
input/manifest checksums. The driver records exact command/resource provenance.
Record selected IDs/reasons, measured elapsed times/bytes/requests and gate outcomes
from the completed run. Never describe the documented SLURM path as a completed
cluster run without scheduler/result evidence.

Retain `context_catalogue.json`, `context_selection.json`, `candidate_pool.tsv`,
`context_distances.tsv`, `context_manifest.tsv`, `combined_metadata.csv`, download
ledger, geography HTML/CSV/TSV/matrix/audit/SVG/PNG, driver provenance and offline
replay audit. Include integrated HTML/supporting archive, native tool logs and
biological gate outcomes; add SLURM job/resource provenance if submitted.

A failed dating gate is valid: retain topology/context and catalogue geography without
claiming a supported dated tree. Full-public, selected-context and focal denominators
remain separate. Composition reflects available sequenced records and surveillance/
submission bias, not prevalence, incidence or transmission. The actual bounded
acquisition/native-analysis result, raw-count check and optional ATB comparison must
be recorded before calling this pilot biologically validated. CI fixture success alone
does not establish those results.

## Checked group and migration comparison

The conflict-aware depth-7 focal prefix has **N=1,767 sample units**, including
Unknown=15 (known-country N=1,752). An independently counted subset is United
States=1,024, Slovenia=213, Italy=181, Oman=29 and Singapore=1; all named-country
counts, including Unknown, sum to N. Counts derive from the full explicitly
QC-passing, accession-deduplicated catalogue before any candidate selection,
not from the selected tree tips.

The limited legacy comparison uses bundled **2025-05** Klebsiella ST147 ATB data:
3,336 already HQ-selected downloadable assemblies with perfect MLST assignment,
2,643 known-country records, and 2,473 dates accepted by the preserved legacy
parser (2,475 under the raw conservative date interpretation). ATB catalogue file
SHA-256 is `43cf66810c4f995d449b90661e4252ffc6980effaad3e12781d5ac1eb494df6b`.
Strong BioSample overlap is 2,873, with 463 ATB-only and 2,848 Pathogenwatch-only
BioSamples. These are unmatched accession sets, not proof of biological absence.
An illustrative seeded country/year-balanced pool of 24 from each provider shares
one BioSample (`SAMN22959219`); this compares pre-download pools, not final
SKA-selected genomes. Provider QC pipelines, upstream sampling and retrieval
times differ. Zero failures in an already HQ-selected ATB catalogue cannot
establish universal upstream ATB QC success; missing raw-failure coverage is
explicitly unavailable. No ATB cgLIN assignment pipeline is implied.

The committed `frozen_catalogue.json.gz` is a permitted compact public freeze with
all 5,647 units and the verified three focal stand-ins. It preserves source UUIDs,
sample units/raw genome counts, normalised country/date/QC and cgLIN annotations,
and source/retrieval hashes, while excluding full raw metadata and assemblies.
Regenerate all figures offline using `scripts/st147_pilot.py --stage figures
--frozen-fixture validation/st147_pathogenwatch/frozen_catalogue.json.gz --output
validation/st147_pathogenwatch/offline_figures`. Its payload hash is checked before
rendering. `test_real_frozen_st147_public_country_and_coverage_hand_counts` verifies
the measured denominators, country subset, assignment coverage and focal overlap.

## Completed bounded preparation

The final live preparation used seed **7**, candidate pool **8**, maximum context
**4**, one nearest candidate per focal, and two threads. All eight candidate
assemblies validated and had SKA comparisons; four were selected. The retained
public focal samples matched their imported cgLIN assignments at depths 5–7, and
were excluded from candidate selection through their accession aliases.

The selected context IDs/reasons are:

| Source genome ID | Country | Selection reason |
|---|---|---|
| aBbUpWyNcgnBayYVXLFaX4 | Italy | nearest to all three focal samples |
| oSavMhBThwcTN6GSL6Rg9q | Italy | balanced background, 2019 |
| 5ifnG8SSwbwRCiwbs5Wahp | Portugal | balanced background, 2009 |
| 9Cg5tSRuj3WekAkohxtNF8 | Nigeria | balanced background, 2016 |

The final ledger materialised **45,436,661 bytes** across eight candidate assemblies:
seven cache hits and one fresh download (two requests, aggregate recorded per-record
time 3.554 seconds). This is not a cold-cache transfer benchmark. Focal acquisition
materialised **17,713,980 bytes**, with six requests and aggregate recorded
per-record time 12.551 seconds. Source-ID mapping, expected provider checksums and
local SHA-256s are retained in the run ledger/manifest.

Audited downstream losses after accession deduplication/focal exclusions were 74
QC-fail/unknown units, 272 undated/invalid units, zero explicit epidemiological
filter losses, 5,364 units outside the bounded pool, zero download/SKA-comparison
failures, and four excluded by the final selection limit. These denominators refer
to selection stages, distinct from the full-public geography denominator.

`frozen_catalogue.json.gz` now includes the four selected source IDs/reasons and
selection seed. Committed `country_composition.csv.gz`, `.tsv.gz` and
`country_composition_matrix.csv.gz` are the final full-public, selected-context and
focal figure source tables. The real-freeze test compares every exported row against
offline regeneration, in addition to the hand-counted subset. The compact fixture's
internal payload SHA-256 is
`9fff151f4d91916ed05169b50b9ecdb3268a53c28142043e2519d0a53d54e3dc`.

## Completed native end-to-end result

The **full** native laptop workflow completed on macOS ARM64 with seven samples,
seven distinct dates, 5,816,827 aligned sites and 5,085,048 retained complete clonal
sites. Raw and clonal lineage-coherence checks flagged **zero of seven** samples.
The root-to-tip temporal-signal gate was **not supported**: observed R²=0, rate
`1.44e-8`, permutation p=1.0 at the unchanged threshold 0.05. All 19 requested
permutations completed. This bounded release plumbing run used 19 permutations;
the reproducible production laptop/SLURM recipe retains 100. The failed gate was
preserved: **zero supported dated trees** were claimed, while topology, geography,
context interpretation, HTML report and supporting archive remained available.

The actual command was:

```bash
chronoclade run tmp/st147-pilot/context/combined_metadata.csv \
  --context-manifest tmp/st147-pilot/context/context_manifest.tsv \
  --output tmp/st147-pilot/results --threads 2 --lineage-jobs 1 \
  --randomisation-jobs 2 --date-randomisations 19 --min-samples 4 --seed 7
```

Recorded tools: SKA 0.5.1; IQ-TREE 3.1.3; ClonalFrameML 1.20; TreeTime 0.12.1.
Evidence comes from actual analysis-log headers and version commands in the same
Pixi environment. Sanitised actual `native_results/` files preserve report/temporal
summaries, raw/clonal coherence rows, screening distances, clonal topology,
execution resources, acquisition checksums and validation ledgers. Machine-local
paths and output pointers are omitted. No dated-tree/node-date result is fabricated
when the temporal gate fails.

The separate small live bulk contract test verified three source-ID/content joins:
each bulk FASTA exactly matched its independently downloaded single FASTA by
SHA-256. Two bulk requests returned 6,253,133 response bytes in 25.77 seconds.
A deliberately partial batch requested one valid source plus invalid numeric ID
`999999999`: exactly one valid FASTA returned and the missing ID remained visible.
All eight candidate assemblies additionally matched their expected provider SHA-1
and sequence-base counts; a subsequent cache-resume validation made **zero network
requests**. Compact public audit files are in `native_results/`.

SLURM submission syntax is validated and a portable staging recipe is supplied.
**No cluster job was submitted** in this pilot. This is a completed native laptop
result plus an exercised offline reproduction path, not a claim of cluster execution.

Representative publication exports are retained in `figures/`: the full-public
prefix-depth-5 percentage overview (page 1), selected-context depth-7 percentages,
and separate focal depth-7 percentages, each as SVG and PNG. They were inspected
at readable size: named countries, full-prefix/scheme labels, N/known/Unknown,
legends, group-size labels, scope and interpretation captions remain legible.
The public overview orders groups by size; every singleton remains in the complete
source tables and paginated offline output. Bounded selected/focal panels show every
named country rather than needlessly collapsing them to Other. Caption/legend
space is reserved separately. The integrated report was also checked in the browser
by the orchestrating agent, including the failed temporal-signal gate.
