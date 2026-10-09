# Output files

The top-level `index.html` links the report for each species/lineage group.
`summary.json` records analysed and skipped groups.

The `logs/` directory contains each external command and a JSON fingerprint of
its direct inputs. ChronoClade uses those fingerprints to decide whether a
completed stage can be reused.

## Main lineage outputs

| File | Content |
| --- | --- |
| `report.html` | Reader-facing analysis report |
| `report.json` | Machine-readable report record |
| `supporting_results.zip` | Portable report evidence bundle |
| `core_alignment.fasta` | SKA2 reference-ordered alignment |
| `lineage_coherence.tsv` | Per-genome raw-distance screen run before IQ-TREE and ClonalFrameML |
| `clonal_lineage_coherence.tsv` | Per-genome recombination-filtered distance screen run before temporal analysis |
| `iqtree.treefile` | Starting maximum-likelihood tree |
| `clonalframeml.labelled_tree.newick` | Recombination-corrected tree |
| `clonalframeml.importation_status.txt` | Recombination intervals by branch |
| `clonalframeml.filtered.fasta` | Alignment after recombination filtering |
| `recombination_intervals.tsv` | Normalised branch-level calls with alignment and reference-record coordinates, including cross-record flags |
| `recombination_genome_profile.csv` | Binned importation depth and filtering values used in the map |
| `recombination_map.svg`, `.png` | Genome-wide recombination and alignment-filtering figure |
| `recombination_summary.json` | Counts of inferred intervals and columns removed or retained |
| `clonal_pairwise_distances.tsv` | Pairwise SNPs, callable sites and comparison classes |
| `clonal_snp_matrix.tsv` | Square clonal SNP matrix |
| `pairwise_callable_sites.tsv` | Square callable-site matrix |
| `public_health_evidence.json` | Scenario rules, groups, neighbours and sensitivity results |

## Temporal evidence

| File | Content |
| --- | --- |
| `clock/root_to_tip_regression.svg` | TreeTime root-to-tip plot |
| `clock/rtt.csv` | Root-to-tip values for every genome |
| `root_to_tip.png` | High-resolution regression plot |
| `temporal_signal.json` | Observed and permuted clock statistics |
| `date_randomisation.csv` | One row per observed or permuted fit |
| `date_randomisation.svg`, `.png` | Randomisation distributions |
| `timetree/timetree.nexus` | Time-scaled tree, when the gate passes |
| `timetree_with_confidence.svg` | Calendar tree with node intervals |
| `timetree.png` | High-resolution dated tree |
| `timetree_confidence.csv` | Clock and root uncertainty summary |
| `node_dates.csv` | Internal-node estimates and 90% intervals |

## Fast-mode evidence

| File | Content |
| --- | --- |
| `phipack_profile.tsv` | PhiPack Profile positions and p-values, with reference-record coordinates |
| `phipack_significant_blocks.tsv` | PHI-positive profile results and their computational block coordinates |
| `phipack_summary.json` | Parameters, tested blocks, detection result and interpretation limits |
| `core_alignment.fasta` | Uncorrected SKA alignment supplied to the fast-mode IQ-TREE run |
| `iqtree_fast.treefile` | Uncorrected screening phylogeny used for root-to-tip analysis |

## Context preparation

`prepare-context` writes the complete same-ST accession list, the balanced
candidate pool, SKA screening distances, the selected context manifest,
combined metadata and a JSON audit of the selection process. Keep these files
with the analysis. A different public context sample can change the apparent
placement of focal isolates.

The Pathogenwatch route also writes the full frozen raw catalogue
`pathogenwatch_catalogue.json`, annotated `context_catalogue.json`,
`download_ledger.json` (successes, failures, checksums, bytes and timing), and
`context_geography/` with count/percentage SVG and PNG figures, source tables,
assignment coverage, deduplication audit and standalone HTML. The public panel
uses the full QC-eligible catalogue; selected-context and focal panels have
separate denominators. These outputs remain usable when dating is unsupported.
Manifests retain source IDs, accessions, cgLIN prefix keys and catalogue hashes.
Keep the referenced frozen catalogue together with the manifest. Lineage reports
copy the annotated catalogue, selection audit and figure tables to supporting
results archives.
