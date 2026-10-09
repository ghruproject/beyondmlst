# ChronoClade

ChronoClade tests for measurable evolution in longitudinal bacterial genome
collections. It starts with a recombination-corrected phylogeny, asks whether
genetic divergence is associated with sampling time, and compares the observed
clock fit with a null distribution made by permuting sample dates.

A positive rate and a date-randomisation p-value at or below the chosen threshold
permit time scaling with TreeTime. Collections that fail this test retain their
clonal phylogeny and distance analysis, but ChronoClade does not report inferred
node dates for them.

## Questions the workflow can address

Within a species and lineage, ChronoClade can help an investigation examine:

- whether isolates sampled months or years apart belong to the same local clade;
- whether local isolates are split into separate parts of the tree;
- whether public genomes change the interpretation of those groups;
- whether the sampling dates support an evolutionary rate and dated tree; and
- how much uncertainty surrounds inferred ancestral dates.

These are population-level questions. The workflow does not reconstruct direct
transmission or count importation events.

## Analysis outline

```text
assemblies and sample dates
          |
          v
whole-genome alignment and maximum-likelihood tree
          |
          v
ClonalFrameML recombination correction
          |
          +---- clonal SNP distances and local/context topology
          |
          v
root-to-tip regression
          |
          v
date-randomisation test
          |
          +---- no temporal signal: stop time scaling
          |
          v
TreeTime dated phylogeny with node-date intervals
```

[Install ChronoClade](installation.md){ .md-button .md-button--primary }
[Follow the ST239 example](worked-example.md){ .md-button }

## Current scope

ChronoClade accepts assembled bacterial genomes. It analyses each
species/lineage combination independently and can use public context selected
from Pathogenwatch, with AllTheBacteria retained as an explicit legacy route. The current implementation uses SKA2, IQ-TREE,
ClonalFrameML and TreeTime in a Pixi environment.

The software is an early working release. Its temporal estimates have been
checked against a published *Staphylococcus aureus* ST239 dataset. Context
selection and public-health interpretation still require review by an analyst
who knows the sampling frame and local epidemiology.
