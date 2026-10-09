# Reading the report

The report starts with what was found, the countries represented and the closest
analysed relatives. Dating comes afterwards: genetic relationships remain useful
when collection dates cannot support reliable estimates of ancestral times.

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

## Fast-screen report

Fast mode remains explicitly limited: recombination detection, a genetic tree
without recombination correction, and a date-permutation screen. It never presents
a dated tree, final clonal-neighbour analysis or an introduction interpretation.
A PHI-positive block signals that the full analysis is needed; it is not itself
a localised recombinant tract.
