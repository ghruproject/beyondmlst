# Public-health questions

## Purpose

`chronoclade` examines longitudinal sampling within a species and ST. Its main
question is whether local isolates collected over months or years are
consistent with one persistent local lineage, several introductions, or a
mixture of both.

The workflow does not infer that one sampled isolate directly descended from
another. A close relationship is evidence of shared ancestry and may include
unsampled intermediates.

## Unit of analysis

Each analysis covers one species and ST, or a narrower genomic lineage when an
ST contains excessive diversity. The user must define the focal population,
such as a hospital, city, region or country. The meaning of local circulation
depends on this definition.

Collection year and month are sufficient inputs. Exact dates should be used
when available. Imprecise dates must retain their stated uncertainty rather
than being converted to arbitrary days.

## Main interpretation

The report assigns one of four evidence summaries:

1. Consistent with a single persistent local lineage.
2. Consistent with multiple introductions.
3. Evidence of local persistence plus additional introductions.
4. Indeterminate because contextual, temporal or genomic resolution is
   insufficient.

These summaries describe consistency with each scenario. They are not counts
of transmission or introduction events.

## Questions the report should answer

### Population structure

- Does the ST contain one local genomic population or several distinct groups?
- How many recombination-corrected local clusters are supported by the data?
- Do isolates from later months or years remain within an earlier local clade?
- Is there evidence of persistence, replacement, disappearance or
  re-emergence?

### Genomic distances

- What are the pairwise clonal SNP distances between isolates?
- How do within-month, within-year and between-year distances compare?
- Are between-cluster distances wider than distances within local clusters?
- How many sites were callable for each comparison?

SNP counts must be calculated after accounting for recombination. The report
should show the observed distributions and callable sites. It must not impose a
universal SNP threshold because appropriate thresholds depend on species,
setting, sampling density and surveillance purpose. A study-specific threshold
may be supplied when it has an external justification.

### Tree topology and context

- Do local isolates form one supported clade or several separated clades?
- Are local groups interspersed with retrospective or public contextual
  genomes?
- Which contextual genomes are closest to each local cluster?
- Are the closest contextual relationships drawn from different locations or
  collection periods?

Contextual genomes from the same ST are needed to distinguish introductions
from local circulation. A local-only dataset can reveal one or several local
clusters, but it cannot establish where those clusters originated.

Candidate context assemblies are fetched from Pathogenwatch for supported
organisms. Selection must include close genomic neighbours of each focal group
and a geography/time-stratified background. A convenience sample based only on
availability or assembly quality is not sufficient. The report must link to a
frozen manifest recording the candidate pool, selection rules, accessions,
metadata snapshot and reasons for inclusion.

### Longitudinal patients and epidemiology

- Do serial isolates from one patient support persistent carriage, gradual
  evolution or reacquisition of a distinct strain?
- Which genetically close isolates come from different patients?
- Do those patients overlap in time and place?
- Could repeated isolates from one heavily sampled patient dominate the
  topology or temporal analysis?

The workflow should provide a sensitivity analysis with one representative
isolate per patient. Ward, hospital and epidemiological overlap can support a
plausible transmission link, but genomics alone cannot establish who infected
whom.

### Time and uncertainty

- Is there sufficient temporal signal to estimate evolutionary dates?
- When did each supported local cluster share a common ancestor?
- Does the inferred ancestor predate the first observed local isolate?
- How wide are the node-date intervals?

The workflow should produce an undated recombination-corrected tree regardless
of temporal signal. It should produce a time-scaled tree only when the
date-randomisation gate passes, and show TreeTime node-date intervals on that
tree.

### Resistance and phenotype

- Is an AMR phenotype or determinant confined to one local cluster?
- Did resistance emerge within a persistent lineage or arrive with a distinct
  introduction?
- Has the resistance profile changed over the sampling period?

Core-genome topology cannot establish plasmid transmission. Plasmid or other
mobile-element questions require a separate analysis, ideally with long-read
or validated plasmid-resolution data.

### Robustness and sampling

- Are the relevant branches supported?
- How much recombination and missing data affect the result?
- Does the conclusion change after removing repeated patient isolates?
- Does it change with a different reference or contextual-genome selection?
- Which additional samples would most reduce the current uncertainty?

## Evidence expected under each scenario

| Evidence | Persistent local lineage | Repeated introductions |
| --- | --- | --- |
| Clonal SNP distances | Small within a local group across the sampling period | Wider distances between distinct local groups |
| Topology | Local isolates remain in one supported clade | Local isolates occupy several separated clades |
| Contextual genomes | Local isolates are closer to one another than to context | Local groups associate with different contextual genomes |
| Time tree | A shared ancestor is compatible with persistence during or before surveillance | Separate ancestral histories predate local detection |
| Location reconstruction | A persistent focal-location state | Several context-to-focal transitions are plausible |

No single row is decisive. Similar introductions can have small distances when
the source population is poorly sampled, and an established local lineage can
accumulate diversity over time.

## Actionable report implementation

The current report now provides:

1. A recombination-filtered pairwise SNP matrix and heatmap, including callable
   sites.
2. Distance distributions within years, between years and between inferred
   clusters.
3. A candidate local-group table listing samples, patients, date span,
   within-group distance and nearest contextual genomes.
4. A deterministic one-isolate-per-patient distance sensitivity view.
5. A cautious scenario summary and explicit evidence ledger supported by
   topology, distances, contextual placement and time estimates where
   available.

The next report iteration still needs a metadata-coloured phylogeny, propagated
branch support, a patient/location timeline and a full topology sensitivity run
with one isolate per patient.

The temporal-signal result remains part of the report, but it is not the main
public-health conclusion.
