# Greece CCRE: Pathogenwatch and REP terminology validation

This fixture uses the existing Greek CCRE focal sets: **36 ST39 and 11 ST147
isolates**, all collected in 2019. It tests public Pathogenwatch identity
resolution, profile coverage and the REP-inspired report wording. The four-query
subsets exercise assembly acquisition and the corrected genomic workflow.
All 47 focal queries have valid fast profiles. Both nominated four-query
subsets completed the native workflow with **100 successful date permutations**
each and withheld time trees because their temporal gates were unsupported.
These are bounded validation datasets, not assessments of each entire lineage.

## Portable inputs and identity evidence

| File | Purpose |
| --- | --- |
| `ST39_queries_public_all.csv` | All 36 ST39 focal isolates |
| `ST147_queries_public_all.csv` | All 11 ST147 focal isolates |
| `ST39_queries_public_subset4.csv` | Four ST39 queries selected using dates and hospitals |
| `ST147_queries_public_subset4.csv` | Four ST147 queries selected using dates and hospitals |
| `exact_accession_mappings.csv` | Exact focal BioSample/run matches and both saved collection and public genome identities |
| `input_audit.json` | Counts, metadata/assembly policies, subset selection and query-file hashes |
| `source_hashes.json` | SHA-256 hashes of the original sources and local input-audit evidence |
| `validation_status.json` | Measured fast coverage, export failures and final native-stage outcomes |
| `ST39_report_output_audit.json` | Existing local HTML targets/anchors and advertised report JSON outputs |
| `ST39_temporal_signal.json` | Actual ST39 root-to-tip/date-permutation results |
| `ST147_report_output_audit.json` | Equivalent ST147 local HTML and JSON output checks |
| `ST147_temporal_signal.json` | Actual ST147 root-to-tip/date-permutation results |
| `ST39_fast_report_output_audit.json` / `ST147_fast_report_output_audit.json` | Final fast-report target checks and observed rendered-report hashes |
| `cgmlst_identity_failure_evidence.json` | Requested/returned public genome IDs and hashes for the rejected export and serial recovery |
| `st147_cgmlst_serial_retry.audit.json` | Exact identity validation of the successful serial focal-export replay |

Every focal BioSample resolves to an ENA run and a saved CCRE collection record.
The replacement public record was verified against **both** BioSample and run
accession; no name-only or fuzzy identity joins were used. Collection and public
genome UUIDs and numeric IDs are retained together so the substitution remains
auditable. Original collection records were omitted by the configured analysis
capability endpoint during this validation.

All query files omit assembly paths. The workflow acquires the corresponding
public Pathogenwatch query assembly. This is a same-isolate alternative assembly
and may differ from the original collection or AllTheBacteria assembly; matching
accessions alone do not establish byte-identical assemblies. Large raw
catalogues, FASTAs, credentials and local assembly paths are excluded from this
fixture. ST39 mappings retain individual raw-detail hashes; ST147 retains the
hash of the complete verified detail response file in `input_audit.json`.

Full publication collection dates, hospital codes and NUTS-2 codes override
year-only public-record dates. `region` is blank because the source supplies
codes, not region display names. `country` and `location` are Greece. No region
names were inferred.

The subset algorithm starts with the earliest collection date, then maximises
minimum date separation, preferring hospitals not already selected. Ties use
date, hospital and BioSample. It uses no genetic group, distance, branch or
clock-test result.

The historical supplied ST147 typing export has 7,804 UUID records and 5,620
deduplicated sample rows. It contains metadata and cgLIN results, **not FASTA
assemblies**. Its source hashes are retained without copying the full export.
The scheme is `scgMLST629_S`; its historical version and export date are unknown.
These historical counts describe a separate pool from the bounded live context.

## Validation snapshot: 9 October 2026

The actual `fast` runs completed for all **47 focal inputs**, with a limit of 80
public context records per lineage and seed `20260927`. Completion preserves
coverage failures rather than implying every focal genome was compared.

| Lineage | Focal profiles available | Context profiles available | Comparable pairs | Fast status |
| --- | ---: | ---: | ---: | --- |
| ST39 | 36 / 36 | 80 / 80 | 6,670 | Completed |
| ST147 | 11 / 11 | 80 / 80 | 4,095 | Completed after export retry |

The first ST147 focal cgMLST export was rejected because it contained an
unrequested genome or an empty locus. That attempt preserved all 11 focal
queries in coverage and exclusions and analysed only 80 context profiles. A
serial export subsequently supplied valid profiles for all 11 focal queries;
the successful retry is shown above. `validation_status.json` retains the
rejected attempt separately. Live export jobs were advertised by the server,
but their typing database version and fingerprint were unavailable. The
rejected ST147 response was byte-identical to the successful ST39 export and
contained the ST39 identities. A transient upstream selection/cache collision
is an inference from that evidence, not a diagnosis of the remote server.

The native ST147 subset acquired four query and six context assemblies, then
failed the lineage-coherence safeguard. Context `PW_73f6h9aWpgxZT9XSRMQNzA` was
flagged; the corrected genomic workflow stopped. The reduced four-context run
subsequently completed all stages: three nearest-profile records and one
genetic-group/time/region representative, selected before clock testing. This
does not validate the default 50-context budget.

ST39 completed `fast`, `full` and `finish` with four focal and six public context
assemblies. Its observed root-to-tip rate was `-1.448e-08`, R² was `0`, and all
100 requested date permutations succeeded with empirical p = `1.0`. ST147's
four-query/four-context run had observed rate `1.474e-06`, R² `0.26`, and 100/100
successful permutations with empirical p = `42/101` (about `0.416`). Neither
subset supported time scaling at the fixed `0.05` threshold, so **no time tree
was produced**. ST39's negative observed rate independently fails the
positive-rate requirement. These results apply to the selected datasets and
do not establish a lineage-wide absence of temporal signal.

Earlier 10-permutation smoke tests also withheld time trees. With 10
permutations the minimum empirical p is `1/11` (about `0.091`), which cannot pass
the `0.05` gate. Those smoke outcomes are retained separately in
`validation_status.json`; the 100-permutation results above supersede them.
The CLI default is 100 permutations. The supported-clock/time-tree branch was
covered by unit tests, not by these native datasets.

The final native validation settings were:

| Setting | ST39 | ST147 final subset |
| --- | ---: | ---: |
| Focal queries | 4 | 4 |
| `--context-size` | 6 | 4 |
| `--nearest-per-query` | 2 | 1 |
| `--profile-limit` | 80 | 80 |
| `--profile-bootstraps` | 10 | 10 |
| `--date-randomisations` | 100 | 100 |
| `--min-samples` | 4 | 4 |
| `--threads` / `--lineage-jobs` / `--randomisation-jobs` | 2 / 1 / 2 | 2 / 1 / 2 |
| `--seed` | 20260927 | 20260927 |

Both used `--mode finish`, `--date-randomisation-method root-to-tip` and
`--temporal-p-value 0.05`.

All local targets and
anchors in 25 HTML reports/fragments and every advertised output in five report
JSONs per lineage exist, including `nearest_neighbours.tsv`.

The ST39 country view counts four focal and five selected public sample units
after QC and accession deduplication within each cohort. Its nine sample-unit
entries across the separately deduplicated cohorts differ from the ten genome
assemblies in the tree; the caption now states this denominator explicitly.
Overlapping sample units can contribute once to each cohort. Regressions check
that duplicate selected records and focal/public overlaps do not become an
incorrect claim about tree size or unique sample counts. See
`validation_status.json` for the machine-readable outcomes. Report hashes were
refreshed after the final style and caption regeneration. Any later rerun or
rendering change requires another audit and new hashes.

The CSV inputs are portable, but exact replay requires the saved live catalogue,
typing exports and assemblies from the validation run, which are retained
outside Git. A fresh public-data run can select a different pool. Source hashes
identify the validated snapshots; they do not make those excluded files
available in this fixture.

## REP-inspired interpretation

Both generated fast reports use **Recurrence and persistence in sampled genomes
(REP-inspired)** and link to [CDC's Reoccurring, Emerging, and Persisting
framework](https://www.cdc.gov/foodborne-outbreaks/php/rep-strains/index.html).
Each marks reoccurring, emerging and persisting epidemiological status as
**not assessed**, with a reason. The report assigns no official CDC REP
designation, including to Klebsiella.

Genetic group stability under locus resampling, observation across collection
years and concentration in time/place remain separate. Selected genome counts
and dated-tip spans cannot establish repeated outbreaks with quiet periods,
increasing illness or its potential, or consistent illness over time. The
existing `temporal_persistence` JSON key remains a descriptive collection-year
summary for compatibility.

Report regression tests cover the unassessed framework, rejection of a
multi-year REP claim, one-year/missing-year limitations and correct rendering
of nearest-relative comparison-group identities. All 17 geography tests pass,
including separate-cohort denominator and overlap regressions. The final
repository suite passed all 360 tests and global lint; documentation and
whitespace checks also passed.
