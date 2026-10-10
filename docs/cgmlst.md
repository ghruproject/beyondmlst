# Independent cgMLST analysis

`cgmlst` reads a validated prepared bundle. It partitions inputs at an explicit
cgLIN depth (5 by default; 6 or 7 are options), includes every matching context
record already in that bundle, and writes one existing-style report per block.
Species, supplied MLST ST and lineage scheme/version/database scope separate
blocks. Declared MLST schemes remain separate. A missing MLST scheme can enter
the sole declared comparison namespace for its species/ST while its metadata
remains unknown; multiple declared schemes make that missing namespace ambiguous
and excluded. It does not automatically narrow to a deeper prefix.

```bash
chronoclade prepare frozen-records.json --catalogues locus-catalogues.json --out prepared
chronoclade cgmlst prepared/dataset.json --out cgmlst --lin-level 5 \
  --subsamples 5 --context-size 50
```

For E. coli with existing compatible HierCC assignments, supply an explicit level:

```bash
chronoclade cgmlst prepared/dataset.json --out ecoli-cgmlst --hiercc-level HC10
```

The command does not yet discover or retrieve fresh public context. Prepare the
input and context records together using existing frozen exports. The partition
index and `partitions.json` distinguish exact included IDs, unmatched context,
unresolved inputs and counts available at LIN levels 5, 6 and 7. These counts are
for the supplied bundle, not the complete public collection. Missing ST or lineage
evidence is reported; a neighbour's label is never copied to an input genome.

Every comparable profile remains in the NJ/ordination analysis and downloadable
tree. Selection does not prune this guide tree. The approved location network and
existing per-block report presentation are reused. Ordination metadata toggles,
large-tree interactive display, improved genetic-group calibration and the
large-dataset backend remain planned work. The current 1,500-record per-block
limit remains enforced; removing it would not make this a 10,000-genome engine.

## Representative alternatives

`--subsamples N` requests alternative representative sets from each block.
All input genomes, explicit pins (`--include ID`) and declared nearest neighbours,
including ties at the requested boundary, are mandatory. The context-size budget
excludes input genomes. Mandatory contexts can exceed it: `budget_overrun`
records the excess rather than dropping required samples.

Nonmandatory representatives vary across reproducible seeds while retaining the
selection strata and quota policy. Every selection records IDs, per-ID reasons,
coverage, missing metadata, source distance units and checksums. Pairwise overlap
and cumulative unique coverage are retained in the ensemble. A block with too few
choices can produce repeated sets; those are explicitly labelled as duplicates,
not presented as independent evidence of convergence.

Each block writes `selection_bundle/ensemble.json` and
`selection_bundle/selections/selection-001.json` (and subsequent alternatives).
The bundle copies its distance evidence and refers to its prepared dataset and
partition manifest using relative paths and hashes. Move the common parent of
these bundles together. The public loaders reject changed sources or members
that do not belong to the exact source partition.

These selection manifests are ready for the planned independent `tree` command;
that command and independent `time`/ensemble comparison are not yet implemented.
Between-selection variation will measure subsampling sensitivity, not a formal
population confidence interval.

## Output ownership

`cgmlst.json` is published after all requested blocks finish. It references each
report and ensemble, the prepared dataset, partition evidence and exact analysis
settings. Existing nonempty output directories are rejected. A failed stage does
not publish a completed manifest. Automatic resume and shared local/SLURM job
execution remain future work; use a new output directory for a revised run.
