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
normalisation is
`f8142089ed06ed4e3dc15fd9816353e693bd26920d98881e4b2883b2362759d4`.
Raw responses are unchanged; re-normalisation removed all 39 apparent country
conflicts because those values differed only by aliases. Use the saved run
manifest/provenance to verify the actual input file hashes.

The full annotated **QC-passing** public catalogue deduplicates to 5,647 sample units.
Measured coverage and group sizes are:

| Full-prefix depth | Assigned units | Unresolved units | Groups | Singleton groups | Largest group |
|---|---:|---:|---:|---:|---:|
| 5 | 5,644 | 3 | 111 | 83 | 4,430 |
| 6 | 5,632 | 15 | 256 | 184 | 2,589 |
| 7 | 5,632 | 15 | 661 | 462 | 1,793 |

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
| k6RA17dLahES7gXCbfqgPv | SAMN21163368 | 2013 | Singapore |
| h1KCajmAVsy9WoF7nmxtTn | SAMN11853676 | 2015 | Oman |
| 7LyY3MxHcNrEftnq7JU8wo | SAMN15868620 | 2017 | United States |

They were selected from QC-passing dated public entries with valid shared depth-7
annotation, different years and strong accession links. This is a public demonstration,
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
