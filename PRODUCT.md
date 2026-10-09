# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Public-health microbiologists and bioinformaticians reviewing bacterial genomes
from longitudinal sampling within a species and lineage such as an MLST sequence
type. They need to decide whether a dated analysis is defensible and whether the
genomic pattern is consistent with local persistence, repeated introductions, or
an unresolved mixture.

## Product Purpose

ChronoClade turns assembled genomes, collection dates, and optional contextual
genomes into a reproducible recombination-aware analysis. Success means a reader
can follow every inferential gate, inspect the supporting evidence, download the
underlying data, and understand what the analysis does and does not establish.

## Positioning

The workflow makes temporal signal a visible prerequisite for time scaling and
combines that decision with recombination-filtered distances, topology, sampling
time, and public context. It does not treat an ST label, a dated tree, or a single
SNP threshold as proof of transmission or introduction.

## Operating Context

The software is a Python command-line workflow managed through Pixi. It produces
a self-contained HTML report plus machine-readable and publication-ready files
for each lineage. Reports may be read offline, printed, shared with collaborators,
or used as supporting material in public-health review.

## Capabilities and Constraints

- Input begins with assembled bacterial genomes and collection dates; read QC and
  assembly are upstream responsibilities.
- Root-to-tip regression is an exploratory first-stage diagnostic, followed by a
  date-randomisation gate before TreeTime time scaling.
- TreeTime is a maximum-likelihood method. Confidence evidence includes clock
  fit, rate uncertainty, node-date intervals, permutation results, and sensitivity
  to sampling and model choices.
- Public context is retrieved from Pathogenwatch, currently for *K. pneumoniae*.
  Selection remains explicit because sampling changes introduction-oriented interpretation.
- Automated public-health scenarios are provisional and require epidemiological
  review. The workflow does not establish direct transmission or exact numbers of
  introduction events.

## Evidence on Hand

- A real public ST239 validation dataset and completed results are stored under
  `validation/st239_baines2015/`.
- The current code produces TreeTime root-to-tip, randomisation, dated-tree, and
  recombination-filtered distance outputs, but the reader-facing report does not
  yet explain or expose them adequately.
- No validated local *Klebsiella pneumoniae* or *Escherichia coli* longitudinal
  case study is currently included.

## Product Principles

1. Explain the scientific question before presenting its statistic.
2. Make every proceed/stop gate explicit and preserve undated results when dating
   is unsupported.
3. Put readable figures and downloadable source data beside each conclusion.
4. State uncertainty and method limitations in plain language.
5. Keep interpretation useful for public health without overclaiming transmission.

## Accessibility & Inclusion

Reports must remain readable on desktop and mobile, printable without losing the
analysis order, navigable with semantic headings and links, and understandable
without relying on colour alone.
