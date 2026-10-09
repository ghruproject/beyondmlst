# Typing new query assemblies

Context discovery starts with the species and seven-locus MLST ST, because these
can be queried in Pathogenwatch. The resulting same-ST catalogue is the possible
comparison pool. Compatible cgMLST profiles and genetic-group assignments then
prioritise a bounded pool, with 25% of its budget reserved for country/year-balanced
background. SKA screens the resulting genomes before the final corrected tree is
built. A code alone is not a substitute for the final genetic comparison.

## Install the tools

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first, then run:

```bash
pixi run chronoclade setup-typing --species klebsiella --species ecoli
```

This installs pinned releases of Pathogenwatch's own tools in an isolated local
directory. `plincer` and `hclink` use separate uv environments. The cgMLST caller
is a Node program, so setup also installs a checksum-verified Node runtime and
its native dependencies. BLAST is provided by ChronoClade's Pixi environment.
The tools do not upload query assemblies.

Software installation does not install reference data. The setup manifest records
which databases are missing, and its query configuration cannot be used until
those resources are ready. The default directory is
`~/.local/share/chronoclade/typing`; use `--destination` for a shared cluster path.

To download and index publicly available cgMLST alleles as well:

```bash
pixi run chronoclade setup-typing --species klebsiella --species ecoli --prepare-db
```

Database downloads and indexing can be substantial. Without database credentials,
the upstream downloader restricts PubMLST/Pasteur data to publicly redistributable
records through 31 December 2024. Current LIN and HierCC reference builds need
their own upstream credentials; a Pathogenwatch API token does not grant this
access. Supply protected files when building the references:

```bash
pixi run chronoclade setup-typing --species klebsiella --prepare-db \
  --pasteur-secrets /protected/pasteur-secrets.json
pixi run chronoclade setup-typing --species ecoli --prepare-db \
  --enterobase-key-file /protected/ecoli-api.key
```

The installed tools and completed databases can be reused without downloading
them for each query. Install the native tools separately on each operating system;
the completed reference files can be shared with a SLURM cluster. Credentials are not included in the
query configuration or analysis outputs.

## Call and refine

Use the query configuration written by setup, once its readiness checks pass:

```bash
pixi run chronoclade type-queries focal.csv \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --output query_typing

pixi run chronoclade prepare-context focal.csv --scheme klebsiella \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --cglin-depth 7 --output context
```

`--typing-config` types both the focal assemblies and the bounded downloaded
public pool with the same frozen databases. This supports new assemblies whose
identifiers are absent from Pathogenwatch. After download, lineage/profile
refinement selects at most `--max-context` genomes for SKA, preserving its
background budget. It does not type every genome in the full public database.

Alternatively, use `--query-typing query_typing/query_typing.json` to import
existing results whose sample identifiers, species and assembly SHA256 hashes
match the focal files. `--public-typing` imports a frozen JSON list or an
`assignments` envelope keyed by exact Pathogenwatch `source_genome_id`. It can
contain cgMLST profiles, cgLIN assignments or HierCC assignments. Importing a
query file alone does not type public assemblies; only compatible public
annotations can refine the pre-download pool in that mode.

## Klebsiella LIN

The [Pathogenwatch cgMLST caller](https://github.com/pathogenwatch-oss/mlst) uses
the `klebsiella_1` analysis scheme. Its output feeds
[plincer](https://github.com/pathogenwatch-oss/klebsiella-lincodes), which links to
the Pasteur `scgMLST629_S` LIN nomenclature. These are distinct names for the
caller and the lineage system. Ambiguous loci are excluded from assignment.
Partial or provisional LIN prefixes remain partial; tied reference profiles
cannot supply a more specific code than their common prefix.

Comparisons use the entire prefix at `--cglin-depth`, with a compatible scheme
and database version. Unversioned codes may compare within one frozen export,
identified as export provenance rather than a database release. Independent
unversioned exports are not silently combined.

## E. coli HierCC

```bash
pixi run chronoclade prepare-context ecoli-focal.csv --scheme ecoli \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --hiercc-level HC1100 --output ecoli-context
```

The `ecoli` discovery shortcut uses Pathogenwatch's `mlst` search field; `mlst2`
is an explicit alternative field. Always check that the chosen ST refers to the
same seven-locus scheme as your input. The `ecoli_1` cgMLST caller pairs with
[hclink](https://github.com/pathogenwatch-oss/hclink) and the EnteroBase
Escherichia/Shigella cgMLST reference. Its HierCC assignments retain every
returned level. `--hiercc-level` selects one level for context prioritisation;
missing levels remain unresolved. HClink references must be at least as recent
as the caller database. The canonical locus order is checked before assignment.

HierCC identifiers describe clusters at named thresholds, not LIN prefixes.
Do not compare HC1100 identifiers with HC10 identifiers or with Klebsiella codes.
The current country-group atlas uses Klebsiella cgLIN; E. coli HierCC currently
supports query assignment and comparison-pool refinement, rather than a HierCC
country-group atlas. Country counts and the location network are organism-neutral.

For cluster execution, install in a cluster-visible directory and run the same
`prepare-context` command inside a batch job. Afterwards, the combined metadata
and context manifest feed ChronoClade's existing SLURM analysis workflow. A live
cluster run is separate from local native validation.

## Evidence and limits

Typing outputs preserve the raw cgMLST and assignment JSON, assembly hash, tool
versions, database hashes, ordered loci, novel allele hashes, missing calls and
assignment status. cgMLST distances use shared callable loci and require at
least 90% overlap across the scheme; they are not SNP distances. Duplicate or
ambiguous calls cannot establish similarity. Conflicting fingerprints, schemes
or versions block comparison.

`context_selection.json` records typing-priority counts, background counts,
per-query comparisons and reasons for same-ST fallback. Unsupported, missing or
unresolved typing never becomes a made-up group. The candidate manifest retains
`pool_selection_reason` separately from the final SKA selection reason. This is
context within the analysed pool, not an exhaustive nearest-neighbour search.
