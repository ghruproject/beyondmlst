# ChronoClade

ChronoClade compares a query collection with public genomic context, then
builds deeper genomic and temporal evidence in stages. The default `fast` mode
uses public typing profiles and metadata; it does not download context
assemblies or run a molecular clock. Use `full` for corrected assembly-based
genomic analysis and `finish` when dates should be assessed.

## Questions the workflow can address

Within a species and lineage, ChronoClade can help an investigation examine:

- which public profiles are closest to the query samples, and how much profile
  coverage supports that comparison;
- which descriptive genomic groups occur in the available profiles;
- whether a group recurs across years (persistence) and, separately, whether
  observations cluster in a particular time-and-place cell (concentration);
- whether corrected genomic evidence changes the interpretation; and
- whether the sampling dates support an evolutionary rate and dated tree.

These are population-level questions. Profile proximity, groups and location
models do not establish direct transmission or count importation events.

## Stages

```text
fast (default): public typing profiles + metadata
  coverage, country/region summaries, PCoA, exploratory NJ, groups,
  persistence across years, time/place concentration, profile neighbours
          |
          v
full: selected context assemblies + all focal assemblies
  recombination-corrected tree, SNP neighbours, country/network evidence
          |
          v
finish: temporal tests and gated time tree
  only when dates and temporal evidence support assessment
```

The `fast` and `full` stages run without a clock. Profile-stage root-to-tip
plots, when available, are exploratory and do not gate or estimate dates. The
location network shows modelled, potentially root-dependent state changes; it
is not evidence of transmission.

[Install ChronoClade](installation.md){ .md-button .md-button--primary }
[Follow the ST239 example](worked-example.md){ .md-button }
[Read the staged workflow guide](staged-workflow.md)

## Current scope

Profile availability and public metadata depend on the configured provider and
account access. The fast report preserves the full catalogue and analysis
subset denominators and reports missing profiles explicitly. `full` and
`finish` require accessible assemblies for the focal samples and selected
context. `finish` is the only public run mode that assesses temporal signal;
it writes a dated tree only when the temporal gate passes. See the
[staged-workflow guide](staged-workflow.md) for current inputs and limitations.

The software is an early working release. Context selection and
public-health interpretation still require review by an analyst who knows the
sampling frame and local epidemiology.
