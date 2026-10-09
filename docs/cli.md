# Command line

## Native query typing

`setup-typing` installs pinned Pathogenwatch tools with uv-managed Python
environments and a native Node runtime. `--species klebsiella` and
`--species ecoli` can be combined. `--prepare-db` downloads/indexes public
cgMLST data; protected `--pasteur-secrets` and `--enterobase-key-file` enable
the corresponding lineage reference builds. Software-only installation writes
a configuration with `ready: false` until the databases are available.

`type-queries METADATA --typing-config CONFIG --output DIRECTORY` runs cgMLST
and the organism's LIN/HierCC assigner locally. `prepare-context` can use the same
`--typing-config` to type focal and bounded public assemblies, or import
`--query-typing` and `--public-typing` files. `--cglin-depth` and
`--hiercc-level` select compatible groups within the initial same-ST pool.
See [query typing](query-typing.md) for examples and database requirements.

## `chronoclade preflight`

Reports whether the native genome-analysis and temporal-analysis tools are
available in the active environment. Profile-first `fast` mode does not invoke
these assembly/tree tools.

```bash
pixi run chronoclade preflight
```

## `chronoclade validate`

Checks metadata, assembly paths, dates, lineage sizes and reference selection.

```bash
pixi run chronoclade validate metadata.csv --min-samples 10
```

## `chronoclade run`

Runs the cumulative profile-first workflow and writes an index for completed
stages. The default `fast` mode stops after profile analysis.

```bash
pixi run chronoclade run metadata.csv [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--output`, `-o` | `chronoclade_results` | Output directory |
| `--threads`, `-t` | `4` | Total CPU budget |
| `--lineage-jobs` | `2` | Maximum concurrent lineages |
| `--randomisation-jobs` | `4` | Maximum TreeTime permutations per lineage |
| `--mode` | `fast` | `fast` profile analysis, `full` corrected genomic analysis, or `finish` temporal analysis |
| `--collection` | none | Pathogenwatch collection UUID or full URL as query input |
| `--accessions` | none | Query accession list or CSV |
| `--species` | none | Declared species for new assemblies or accession lookup |
| `--catalogue` | none | Frozen public metadata catalogue |
| `--public-typing` | none | Frozen public cgMLST and group assignments |
| `--query-typing` | none | Verified query typing JSON |
| `--typing-config` | none | Native caller and reference configuration for query assemblies |
| `--cglin-export` | none | Frozen public cgLIN assignments |
| `--profile-limit` | `500` | Maximum public records in profile analysis; query records are additional |
| `--context-size` | `50` | Maximum selected context assemblies per lineage; query assemblies are additional |
| `--nearest-per-query` | `3` | Nearest profile relatives considered per query during context selection |
| `--profile-bootstraps` | `10` | Requested locus-bootstrap replicates for descriptive group support |
| `--group-distance` | `0.02` | Exploratory complete-linkage cgMLST mismatch fraction |
| `--date-randomisations` | `100` | Number of tip-date permutations in `finish` |
| `--date-randomisation-method` | `root-to-tip` | `root-to-tip` screen or `full-tree` TreeTime refits in `finish` |
| `--temporal-p-value` | `0.05` | Temporal gate threshold |
| `--min-samples` | `10` | Minimum genomes per lineage for whole-genome workflow |
| `--seed` | `20260818` | Selection and analysis seed |
| `--context-manifest` | none | Frozen manifest from `prepare-context` |
| `--force` | false | Rerun completed stages |
| `--dry-run` | false | Validate and print the plan only |

`--mode fast` compares the query set with the available public typing context.
It does not download context assemblies or build a corrected tree. Its report
shows query/context profile coverage, exclusions, country figures, PCoA and
neighbour-joining views where possible, descriptive genetic groups, separate
year-persistence and time/place-concentration summaries, an exploratory
allele-distance root-to-tip diagnostic, and nearest profile relatives.

`--mode full` selects up to `--context-size` public context assemblies per
lineage, retains all query genomes, and runs the corrected genomic analysis.
It stops before date tests and time scaling. `--mode finish` runs the same
corrected analysis and then performs the temporal tests; TreeTime creates a
dated tree only when the configured temporal gate passes. If there are fewer
than three usable distinct collection dates, `finish` retains the corrected
report and records that temporal analysis was not assessed.

Supply one query source: a metadata CSV, `--collection`, or `--accessions`.
`--dry-run` writes a plan without fetching typing exports, downloading
assemblies or running analyses. `--public-typing` and `--query-typing` import
local frozen files; they do not grant access to protected production databases.

## `chronoclade prepare-context`

Builds a reproducible public context set for one species and ST.

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme SCHEME \
  --st ST \
  --output DIRECTORY [OPTIONS]
```

Use `pixi run chronoclade prepare-context --help` for metadata filters and
screening controls. Start with `--dry-run` to inspect the candidate pool before
downloading assemblies.

Pathogenwatch is the built-in provider, currently supporting *K. pneumoniae*
with `--scheme klebsiella`. Other organisms and schemes return an explicit
unsupported-provider error. `--catalogue FILE` imports a verified frozen JSON
catalogue without live metadata access. Repeated runs reuse the
frozen metadata; `--refresh-catalogue` explicitly fetches a new snapshot.
`--cglin-export FILE` imports a validated CSV/TSV/JSON export and
`--focal-crosswalk FILE` imports sample IDs and strongly verified public IDs.
Without focal assignments, the report states that limitation.

Selection is recomputed when called; downloads resume from checksum-verified
source-ID caches.

Country figures can be regenerated offline from `context_catalogue.json`
without running phylogeny. See [context geography](context-geography.md).
