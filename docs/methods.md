# Methods

## Whole-genome phylogeny

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

## Root-to-tip analysis

TreeTime reroots the recombination-corrected phylogeny by least squares and fits
root-to-tip distance against sampling date. ChronoClade disables automatic
clock-outlier removal at this gate so that samples are not discarded because
they weaken the temporal fit.

Root-to-tip regression is useful for visual inspection. It is not treated as a
formal temporal-signal test because population structure can associate genetic
distance with sampling date even when the dataset cannot estimate evolutionary
times reliably.[^murray]

## Date-randomisation test

The default `root-to-tip` method permutes sampling dates among tips while
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

## Fast screen

Fast mode uses PhiPack's Profile program with Parsnp's 100-site PHI window,
100-site step and p < 0.01 threshold. Because SKA writes reference records
consecutively, ChronoClade divides each record into bounded 250 kb analysis
blocks; a PHI calculation cannot cross a contig join. Ambiguous and missing
calls remain inside their reference segment and are handled by PhiPack.

A PHI-positive result means that incompatibility was detected somewhere within
a computational block. It does not identify a tract, and ChronoClade does not
mask sites from the fast alignment. The per-block p < 0.01 threshold is
unadjusted: this favours sensitivity for triage, so a positive result is an
escalation signal rather than a tract call or a final recombination inference.
The fixed blocks bound runtime but are not biological segments or Parsnp locally
collinear blocks. IQ-TREE therefore builds the fast screening tree from the
uncorrected SKA alignment. Fast-mode output is suitable for triage, not final
dating or public-health interpretation.

## Time scaling

Passing datasets are analysed with TreeTime using marginal time inference,
covariation and 90% confidence intervals. The output tree has branch positions
on a calendar-time axis. ChronoClade exports the root estimate, clock-rate
uncertainty and each internal node's date interval.

The node intervals quantify uncertainty within the fitted TreeTime model. They
do not account for uncertainty caused by incomplete sampling, metadata error or
an inappropriate clock model, so those limitations must remain visible in the
interpretation.

## Location states and public context

TreeTime's mugration model reconstructs the supplied `location` state on the
dated tree, or on the corrected genetic tree when temporal signal is absent.
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
does not run this reconstruction.

The context workflow freezes a full public Pathogenwatch same-ST catalogue with
source IDs, accessions, QC, raw metadata, date precision and retrieval hashes.
It currently supports *K. pneumoniae*; other organisms and schemes return an
explicit unsupported-provider error. Genome availability does not imply cgLIN
support for other organisms.
QC-passing records are deduplicated by accession-supported sample units,
without claiming unique patients or infections. Undated records remain in
country summaries but are excluded from the dated analysis cohort. It balances the candidate pool
across country and year, screens candidates against every focal sample with SKA
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

## References

- Sagulenko P, Puller V, Neher RA. 2018. [TreeTime: Maximum-likelihood phylodynamic analysis](https://doi.org/10.1093/ve/vex042). *Virus Evolution* 4:vex042.
- [StrainHub project](https://github.com/abschneider/StrainHub). ChronoClade's country view is stylistically inspired by this tool but uses TreeTime's location reconstruction.
- Didelot X, Croucher NJ, Bentley SD, Harris SR, Wilson DJ. 2018. [Bayesian inference of ancestral dates on bacterial phylogenetic trees](https://doi.org/10.1093/nar/gky783). *Nucleic Acids Research* 46:e134.
- Didelot X, Wilson DJ. 2015. [ClonalFrameML: efficient inference of recombination in whole bacterial genomes](https://doi.org/10.1371/journal.pcbi.1004041). *PLoS Computational Biology* 11:e1004041.
- Bruen TC, Philippe H, Bryant D. 2006. [A simple and robust statistical test for detecting the presence of recombination](https://doi.org/10.1534/genetics.105.048975). *Genetics* 172:2665-2681.
- Murray GGR et al. 2016. [The effect of genetic structure on molecular dating and tests for temporal signal](https://doi.org/10.1111/2041-210X.12466). *Methods in Ecology and Evolution* 7:80-89.

[^murray]: Murray et al. showed that temporal and genetic structure can inflate root-to-tip regressions and date-randomisation tests.
[^duche]: A clustered permutation should be considered when collection date is confounded with a well-supported genetic group.
