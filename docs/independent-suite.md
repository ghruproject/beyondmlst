# Independent analysis modules

ChronoClade saves each stage separately. Choose **cgmlst** or the optional
**esm2** route after preparation, then pass one saved selection to **tree** and
its saved phylogeny to **time**. The commands retain the established report style
and share the country state-change network component.

## Prepare inputs

```bash
pixi run -e prepare chronoclade prepare records.json --input-kind profiles \
  --catalogues catalogues.json --out prepared

pixi run -e prepare chronoclade prepare \
  https://next.pathogen.watch/collections/jX5cwsoyJ1KDquMqssUzAD-ellington-et-al-2015 \
  --input-kind collection --out collection-input
```

The collection becomes the input dataset. Preparation does not discover public
context. It records profile and lineage readiness, exact identity links, missing
metadata and provenance. Exact linked ENA metadata can fill missing host, country,
collection date and isolation source; conflicting values remain inspectable.

Assemblies require an explicit species and configured typing tools **and reference
databases**. Installed callers do not make an unavailable database ready. Existing
valid profiles and assignments can be imported while reference access is pending.
Credentials are read from runtime configuration; do not put them in input files.

## Explore the complete compatible profile block

```bash
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --out cgmlst \
  --lin-level 5 --subsamples 5 --context-size 50
```

Use `--fetch-context` to discover and retrieve compatible public profiles and
metadata first. The initial pool uses species/ST, followed by the explicitly
requested LIN depth (5, 6 or 7) or HierCC level. Context assemblies are not
requested at this stage. The retrieval audit distinguishes discovered,
retrieved and usable records; an observed export is not a complete reference
locus catalogue.

The partition index links each disjoint block report and its selections. Trees
and ordination use the complete comparable profile block. Ordination views switch
between dataset, country and date colours without changing coordinates; marker
shape distinguishes input from public comparison genomes in every view.

Local groups use an explicit complete-linkage distance threshold, with threshold
sensitivity and optional locus-resampling support. These are descriptive genetic
groups, not automatic outbreak or expanding-clone classifications. Nearest
neighbours report raw allele differences, compared loci and ties.

Multiple selections are alternative representative sets from the same block.
Every input, pinned context and required closest relative is retained, including
ties that exceed the requested context budget. Duplicate alternatives are labelled
when the pool offers insufficient choice.

## Build a selected assembly tree and date it

```bash
pixi run -e tree chronoclade tree \
  cgmlst/BLOCK/selection_bundle/selections/selection-001.json --out tree-001

pixi run -e time chronoclade time tree-001/tree.json --out time-001
```

Tree acquires assemblies only for the exact selected IDs. It runs alignment,
phylogeny and recombination adjustment, then saves temporal evidence and its own
report. Failed acquisition does not choose replacements. Time reads the saved
tree, alignment/site-count and assessment evidence without rebuilding the tree.
Unsupported temporal evidence produces an explicit result rather than a dated
tree. Any supported override is explicit and recorded.

Repeat with the other saved selections, then compare:

```bash
pixi run -e time chronoclade time-compare \
  cgmlst/BLOCK/selection_bundle/ensemble.json \
  --run time-001/time.json --run time-002/time.json --out comparison
```

The comparison uses exact shared anchors and target pairs. It distinguishes the
shared-anchor ancestor date from whole-selection roots, which may describe
different ancestors. Between-selection ranges measure sampling sensitivity;
they are not confidence intervals or proof of convergence.

## Optional embeddings

The ESM2 route requires genuine scoped allele sequences, a validated mapping and
an explicitly acquired checkpoint or explicitly permitted model download. Protein
embeddings cannot observe synonymous changes. Distances and date slopes have
embedding units, not substitutions/site or mutation-rate units. The conventional
cgMLST baseline supplies the biological guide tree and country reconstruction.
See [ESM2](esm2.md) and [sequence inputs](esm2-sequence-input.md).

## Local and SLURM jobs

The same stage command can be saved as a local or SLURM job:

```bash
chronoclade job --slurm --account YOUR_ACCOUNT --partition YOUR_PARTITION \
  --cpus 8 --memory-gb 32 --out tree-job tree -- \
  cgmlst/BLOCK/selection_bundle/selections/selection-001.json --out tree-001 --threads 8
```

`job` writes a quoted script and a command/resource manifest. It does not submit.
Set the executable to an installed stage environment's `chronoclade` path when
necessary. Transfer the dataset, selections and their referenced files together,
and configure site resources and credentials independently. Script generation
and local testing do not establish successful execution on a real SLURM cluster.
