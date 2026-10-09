# Methods

## Stages and profile analysis

The public `chronoclade run` modes are cumulative. `fast` (the default) analyses
available cgMLST profiles without downloading public assemblies. `full` selects
context assemblies and runs the recombination-corrected genomic workflow but
stops before clock fitting. `finish` adds the temporal-signal tests and gated
time scaling. The profile-stage root-to-tip diagnostic and the finish-stage
clock test answer different questions; the former never passes the temporal
gate.

Each profile comparison requires compatible species, cgMLST scheme/version,
database provenance where available, and a compatible locus set. Allele
identifiers are categorical; numeric allele IDs are not evolutionary
distances. Pairwise distances are the number of mismatches divided by the
number of jointly called loci. Missing calls are excluded from the numerator
and the report gives the shared-call count and overlap. Comparisons with
incompatible schemes or locus sets are not combined into one cohort.

Some exports list only observed loci and do not declare the full canonical
scheme universe. For such compatible records, ChronoClade uses the union of
observed loci as the overlap denominator and flags that canonical completeness
is unknown. This does not estimate the proportion of the full scheme covered.
Complete-comparability cohorts are deterministic partitions for ordination and
neighbour-joining views; they are not themselves genetic groups.

Genetic groups use deterministic complete linkage at the configured
`--group-distance`, expressed as a mismatch fraction. Locus-bootstrap
co-assignment is descriptive support for those candidate groups. If no valid
bootstrap replicates are available, support is unknown, not evidence against a
group. A singleton has no within-group pair support to estimate.

Genetic group stability under locus resampling and concentration in time/place
are reported separately.

The machine-readable collection-year summary lists the observed years in the group; gaps do not
demonstrate uninterrupted presence. Concentration summarizes year/location
cells among records with both date and place metadata; it does not test for an
outbreak or establish transmission. Both summaries describe the sample
submitted, not population prevalence. The existing `temporal_persistence` JSON
key retains this descriptive collection-year meaning for compatibility.

The profile-stage nearest relative is selected from the available compatible
public profiles and reports shared called loci. Country figures include
profile-available members of a comparable cohort, combine query and public
context records, and exclude records without profiles. They do not split
country counts by origin or provide a separate region-count view in the current
profile summary. The separate query/context coverage counts include records
without profiles.

## Whole-genome phylogeny

The following assembly methods apply to `full` and `finish`. Clock methods
apply only to `finish`.

ChronoClade analyses each `species` and `lineage` group independently. SKA2
builds a split-k-mer index and maps samples to a selected reference to produce a
reference-ordered alignment. Before tree inference, ChronoClade compares each
genome's median raw SNP proportion with the cohort distribution. A genome is
flagged only when it is both at least five times the cohort median and an
extreme robust outlier (robust z-score at least 10). The workflow stops and
writes `lineage_coherence.tsv` when this screen fails, avoiding an expensive
tree and recombination run on an obvious accession, species or lineage error.
After ClonalFrameML, the same relative check is repeated on the filtered
alignment and written to `clonal_lineage_coherence.tsv` before clock analysis.
This relative check is an input safeguard, not a universal bacterial SNP
threshold or a transmission definition.

IQ-TREE estimates the starting maximum-likelihood phylogeny with a GTR+G model.

ClonalFrameML estimates recombination on that tree and reports 1-based, closed
alignment intervals for each affected branch. With `-output_filtered true` and
`-ignore_incomplete_sites true`, it constructs one shared clonal alignment by
removing every column assigned to an import on any branch and every column with
an ambiguous base in any sequence. ChronoClade verifies that these two masks
exactly account for the filtered-alignment length. It exports the original calls,
a normalised interval table, a binned genome profile and SVG/PNG figures.
Reference-record coordinates are recovered from the ordered reference FASTA.
Intervals that cross a record boundary are retained for audit but flagged because
the adjacency was introduced by concatenation rather than by the chromosome.

A branch-level interval is evidence that ClonalFrameML assigned an import to
that lineage of the tree; it does not mean that every sampled genome acquired
the segment. Pairwise clonal SNP counts use A, C, G and T sites callable in both
members of each pair. The report always gives the callable-site count beside the
SNP count. These corrected distances provide a second opportunity to identify
unusual genomes after recombination has been modelled.

## Finish-stage root-to-tip analysis

TreeTime reroots the recombination-corrected phylogeny by least squares and fits
root-to-tip distance against sampling date. ChronoClade disables automatic
clock-outlier removal at this gate so that samples are not discarded because
they weaken the temporal fit.

Root-to-tip regression is useful for visual inspection. It is not treated as a
formal temporal-signal test because population structure can associate genetic
distance with sampling date even when the dataset cannot estimate evolutionary
times reliably.[^murray]

## Date-randomisation test

The `finish` stage's default `root-to-tip` method permutes sampling dates among tips while
retaining the corrected tree and sequence length. TreeTime repeats the same
least-squares clock analysis for each permutation. The empirical one-sided
p-value compares observed R² with the permutation distribution and includes a
one-count correction. R² drives this screen; the accompanying rate distribution
is descriptive.

The optional `full-tree` method reruns the complete TreeTime dating fit for every
permutation and re-estimates the root each time. ChronoClade reports the weak CR1
and strict CR2 criteria used in bacterial tip-date randomisation tests. Its
configured decision uses CR2: the observed approximate 95% rate interval must
overlap none of the corresponding intervals from randomised-date fits. These are
normal approximations from TreeTime's reported rate standard error, not Bayesian
credible intervals.

A dataset passes when:

1. every requested randomisation completes;
2. the observed rate is positive; and
3. the corrected empirical p-value is at or below the configured threshold.

The default is 100 permutations and p <= 0.05. This is a pragmatic gate for the
workflow, not a universal definition of temporal signal. Structured datasets
may require clustered permutations or sensitivity analyses that preserve known
population groups.[^duche]

## Stage boundaries

`fast` does not align assemblies, screen recombination with PhiPack, build a
whole-genome tree or run date permutations. It reports the profile comparisons
and metadata available at that stage. `full` builds the corrected genomic
evidence but does not run the clock methods above. Only `finish` assesses
temporal signal. With fewer than three usable distinct collection dates,
`finish` records that temporal analysis was not assessed and leaves the
corrected report available; it does not create a failed permutation result.

## Time scaling

Passing datasets in `finish` are analysed with TreeTime using marginal time inference,
covariation and 90% confidence intervals. The output tree has branch positions
on a calendar-time axis. ChronoClade exports the root estimate, clock-rate
uncertainty and each internal node's date interval.

The node intervals quantify uncertainty within the fitted TreeTime model. They
do not account for uncertainty caused by incomplete sampling, metadata error or
an inappropriate clock model, so those limitations must remain visible in the
interpretation.

## Location states and public context

TreeTime's mugration model reconstructs the supplied `location` state on the
dated tree when time scaling succeeds, or on the corrected genetic tree when
dating was unsupported or not assessed in `full`.
This reconstruction is exploratory. Its result depends on how locations were
defined and sampled.

The full report's StrainHub-style location network groups sampled tips by the
`location` values in metadata and count, then draws representative changes
between ancestral states on the rooted analysed tree. Labels are used as
supplied: there is no automatic geography hierarchy lookup. Use country labels
for a country-level view; city or hospital labels produce a network at that
level. Unknown or missing location values are excluded rather than passed to
the model as a location state. The view uses the existing TreeTime
maximum-likelihood reconstruction; it does not use StrainHub's parsimony
algorithm. Solid arrows
require unique marginal state assignments at both endpoints with confidence of
at least 0.9. Dashed arrows indicate a tied/ambiguous endpoint or confidence
below 0.9. The lower endpoint confidence is a display rule, not a joint
transition probability. Marginal assignments and changes are conditional on
the rooted tree, sampling and supplied location states. They do not prove
transmission or migration, and they do not estimate national prevalence.
Because reconstruction uses the rooted genetic tree, the network can be shown
after dating is unsupported; it contains no inferred timing. Fast-screen mode
does not run this reconstruction. Profile `fast` has a separate exploratory
maximum-parsimony network based on possible changes across alternate sample
roots; its root fraction is a sensitivity summary, not a probability.

The context workflow freezes a full public Pathogenwatch same-ST catalogue with
source IDs, accessions, QC, raw metadata, date precision and retrieval hashes.
It supports *K. pneumoniae* and *E. coli*; other organisms and schemes return an
explicit unsupported-provider error. Klebsiella uses cgLIN; E. coli uses HierCC.
QC-passing records are deduplicated by accession-supported sample units,
without claiming unique patients or infections. Undated records remain in
country summaries but are excluded from the dated analysis cohort. Same-ST
membership defines the possible pool; compatible cgLIN prefixes, HierCC clusters
and cgMLST distances prioritise candidates, reserving 25% of the bounded pool for
country/year background. Unresolved or incompatible typing falls back explicitly
to balanced same-ST sampling. Native typing can process query assemblies and the
bounded downloaded public pool with the same databases before SKA. It screens
the resulting candidates against every focal sample with SKA
distance, then retains nearby genomes and a stratified background. The manifest
records the bounded search. "Nearest" means nearest within the downloaded pool,
not nearest among all public bacterial genomes.

Country proportions describe available sequenced public records within each
scheme/version-scoped full cgLIN prefix. Unknown country contributes to the
denominator. Incomplete/provisional codes retain their status and missing codes
are assignment coverage, not a single lineage. Prefix depth is not a pairwise
SNP threshold. These figures do not estimate prevalence, incidence, migration,
transmission or country of acquisition. Submission bias, missing annotations and
date uncertainty constrain interpretation. See [cgLIN](cglin.md) and
[geography](context-geography.md).

Native typing uses Pathogenwatch's cgMLST caller and `plincer` or exact `hclink`.
Missing, ambiguous and novel alleles, incomplete assignments and database
provenance remain explicit. HClink is a reference-linking inference rather than
an official new EnteroBase cluster designation. See [query typing](query-typing.md).

## References

- Sagulenko P, Puller V, Neher RA. 2018. [TreeTime: Maximum-likelihood phylodynamic analysis](https://doi.org/10.1093/ve/vex042). *Virus Evolution* 4:vex042.
- [StrainHub project](https://github.com/abschneider/StrainHub). ChronoClade's country view is stylistically inspired by this tool but uses TreeTime's location reconstruction.
- Didelot X, Croucher NJ, Bentley SD, Harris SR, Wilson DJ. 2018. [Bayesian inference of ancestral dates on bacterial phylogenetic trees](https://doi.org/10.1093/nar/gky783). *Nucleic Acids Research* 46:e134.
- Didelot X, Wilson DJ. 2015. [ClonalFrameML: efficient inference of recombination in whole bacterial genomes](https://doi.org/10.1371/journal.pcbi.1004041). *PLoS Computational Biology* 11:e1004041.
- Bruen TC, Philippe H, Bryant D. 2006. [A simple and robust statistical test for detecting the presence of recombination](https://doi.org/10.1534/genetics.105.048975). *Genetics* 172:2665-2681.
- Murray GGR et al. 2016. [The effect of genetic structure on molecular dating and tests for temporal signal](https://doi.org/10.1111/2041-210X.12466). *Methods in Ecology and Evolution* 7:80-89.

[^murray]: Murray et al. showed that temporal and genetic structure can inflate root-to-tip regressions and date-randomisation tests.
[^duche]: A clustered permutation should be considered when collection date is confounded with a well-supported genetic group.
