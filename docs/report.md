# Reading the report

The report starts with what was found, the countries represented and the closest
analysed relatives. Dating comes afterwards: genetic relationships remain useful
when collection dates cannot support reliable estimates of ancestral times.

Genome labels prefer a sequencing-run accession (SRR, ERR or DRR), then a
BioSample accession, then an assembly accession when available. User-supplied
focal sample names are retained. Internal Pathogenwatch IDs remain in the
machine-readable results and label mapping; repeated accessions are disambiguated
so separate genome records remain identifiable. The fast report includes
`sample_labels.csv`, and assembly-stage snapshots retain `sample_labels.json`.

## Start with the summary

The opening summary separates your focal samples (the genomes you want to
investigate) from context genomes (the comparison set). It states whether dating
is supported. An unsupported dating result means the analysis completed but did
not establish enough clock information; it does not mean the software failed.

Public-data demonstrations are labelled explicitly. Their focal genomes are
stand-ins, not a defined outbreak or patient cohort, and their automated
interpretations must not be treated as epidemiological findings.

## Countries and genetic groups

For Pathogenwatch runs with a selection audit, the report shows the matching
record count, eligible sample count, screened pool and final comparison count
before the country figures. Expand the filtering steps to see quality checks,
duplicate removal, focal exclusions, usable-date requirements and metadata
filters. The download pool is sampled across country/year groups with a fixed
seed; closest screened candidates are retained first, with shared neighbours
counted once, and remaining places filled by background comparisons. Pool and
comparison limits are run settings, not fixed biological thresholds.

Country counts first compare focal samples with the context genomes selected for
the tree. A separate figure describes public genomes in the cgLIN groups matching
the focal samples. cgLIN codes describe genetic groups at increasingly specific
levels; they are not SNP-distance cutoffs or transmission clusters.

Keep the denominators separate. A tree may contain seven analysed genomes while
the public catalogue contains thousands of quality-checked, deduplicated sample
units. Country percentages within a public group describe that available group,
not the selected tree participants. Unknown countries remain in the denominator.
The full catalogue atlas and source tables are available separately.

These figures describe submitted sequence records. Surveillance and submission
coverage affect the proportions; country metadata does not establish where an
infection was acquired, migration direction or population prevalence.

## Closest relatives and the genetic tree

The first tree shows genetic relationships before time scaling. Its branch
lengths represent model-estimated substitutions per site, not years. Focal tips
are distinguished from comparison genomes and labelled with country and date.

Per-sample neighbour tables use the final recombination-filtered SNP comparisons.
A SNP is a differing DNA base at a position that can be compared between the two
genomes. The callable-site count states how much sequence was comparable. Pairs
with no comparable sites cannot establish closeness. Exact distance ties remain
visible; the nearest rows are not a unique winner when several genomes tie.

A separate ranking uses patristic distance: the total length of branches linking
two tips in the genetic tree. SNP counts and tree distances can rank relatives
differently because they are different measures. Neither requires a dated tree.
The earlier SKA distances were used to select a screening pool; they are not the
final recombination-filtered relatedness result.

“Closest” always means among the genomes analysed in this run. A bounded context
selection cannot establish the globally nearest public genome. Close genetic
relationships alone do not establish transmission or infection direction.

## Country network

The full report may include a StrainHub-style view of location states on the
rooted analysed tree. It uses the `location` values in the analysis metadata as
supplied; there is no automatic lookup or conversion through a geography
hierarchy. Supply country labels for a country-level network. City, hospital or
other location labels produce a network at that level. Nodes summarise the
locations represented by sampled tips and their sample counts. Arrows summarise
representative ancestral state changes in the tree; they are not proven
transmission links, migration routes, country of acquisition or estimates of
national prevalence. Unknown or missing location values are excluded from the
network rather than treated as a location state.

TreeTime assigns marginal probabilities to ancestral location states. A solid
arrow means both endpoint states have a unique assignment with probability at
least 0.9. Dashed arrows mark an ambiguous assignment or an endpoint below that
threshold. Endpoint confidence is not the joint probability of the transition. Ties between states remain ambiguous. These
probabilities are conditional on the rooted tree, its sampling and the location
states supplied to the model; they do not capture uncertainty in those inputs.
The view uses TreeTime's maximum-likelihood ancestral reconstruction, not
StrainHub's parsimony algorithm. It is available after an unsupported dating
result because it uses the rooted genetic tree and does not require dates.

Fast-screen reports do not reconstruct ancestral location states and do not
include this network.

## Interpretation

The working interpretation combines the corrected topology, SNP comparisons,
sampling span and comparison genomes. Detailed evidence states what is available
and what is missing. Automated labels are provisional: compare them with patient
movement, referral, travel and sampling information.

| Label | Sampled pattern |
| --- | --- |
| Persistent local lineage | One focal-only group spans more than one sampling date |
| Multiple introductions | Focal isolates form at least two separated groups and none spans dates |
| Mixed | Several focal groups are present and at least one spans dates |
| Indeterminate | Focal, longitudinal or contextual evidence is insufficient |

## Can this tree be placed on a calendar?

The report then explains whether genetic divergence increases with collection
time and whether the observed relationship is stronger than results with shuffled
dates. In the default root-to-tip screen, R² describes the fit and the empirical
p-value compares it with the shuffled results. R² alone is not a pass criterion.
Full TreeTime refits use the separately documented rate-interval criterion.

Passing datasets receive a second view: the genetic tree placed on a calendar
axis by TreeTime. This estimates ancestral dates under a clock model; it is not
an independent confirmation of the topology or a transmission tree. Date
intervals are conditional on the model, sampling dates and supplied topology.
Unsupported datasets retain their genetic tree and neighbour analyses, with no
endorsed ancestral calendar dates.

## Methods and downloads

Detailed alignment checks, recombination maps, coordinate tables and clock
statistics are expandable. Raw and corrected distance screens detect grossly
divergent genomes; they do not prove every accession or metadata value is correct.

ClonalFrameML's shared filtered alignment removes a column from every genome when
an import is inferred on any branch. Ambiguous-base filtering is counted
separately. Multi-record reference joins are flagged for review because artificial
adjacency can affect inferred intervals.

Figures are exported as SVG and PNG, numerical values as CSV/TSV, trees as
Newick/Nexus and decisions as JSON. The supporting-results ZIP holds the evidence
actually available for that run.

## Staged reports

The public `fast` report is profile-first: it describes profile coverage,
country/region composition, exploratory PCoA and neighbour-joining views,
descriptive genetic groups, profile neighbours, observations across collection
years and time/place concentration. **Recurrence and persistence in sampled
genomes (REP-inspired)** explains [CDC's REP terminology](https://www.cdc.gov/foodborne-outbreaks/php/rep-strains/index.html)
and marks reoccurring, emerging and persisting status as not assessed, with
reasons. The report assigns no official CDC REP designation, including to
Klebsiella. Genetic group stability under locus resampling, observations across
years and time/place concentration answer separate questions. Selected genome
counts and dated-tip spans cannot establish illness trends, outbreaks or
consistent illness. Profile-stage root-to-tip diagnostics and location networks
are exploratory, not clock gates or transmission evidence.

`full` produces corrected assembly-based genomic evidence but does not assess
dates. `finish` adds date randomisation and a dated tree only when the temporal
evidence gate passes. See [the staged workflow](staged-workflow.md) for the
inputs, outputs and limitations of each mode. Older PhiPack screening files,
where present in preserved results, belong to the legacy per-lineage fast
screen and are not outputs of the current public `fast` mode.
