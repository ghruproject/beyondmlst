# Command line

## `chronoclade preflight`

Reports whether SKA2, IQ-TREE, ClonalFrameML, PhiPack Profile and TreeTime are available in the
active environment.

```bash
pixi run chronoclade preflight
```

## `chronoclade validate`

Checks metadata, assembly paths, dates, lineage sizes and reference selection.

```bash
pixi run chronoclade validate metadata.csv --min-samples 10
```

## `chronoclade run`

Runs all ready lineages and writes the combined index.

```bash
pixi run chronoclade run metadata.csv [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--output`, `-o` | `chronoclade_results` | Output directory |
| `--threads`, `-t` | `4` | Total CPU budget |
| `--lineage-jobs` | `2` | Maximum concurrent lineages |
| `--randomisation-jobs` | `4` | Maximum TreeTime permutations per lineage |
| `--date-randomisations` | `100` | Number of tip-date permutations |
| `--mode` | `full` | `full` analysis or `fast` triage screen |
| `--date-randomisation-method` | `root-to-tip` | `root-to-tip` screen or `full-tree` TreeTime refits |
| `--temporal-p-value` | `0.05` | Temporal gate threshold |
| `--min-samples` | `10` | Minimum genomes per lineage |
| `--seed` | `20260818` | Randomisation seed |
| `--context-manifest` | none | Frozen manifest from `prepare-context` |
| `--force` | false | Rerun completed stages |
| `--dry-run` | false | Validate and print the plan only |

`--mode fast` always uses `root-to-tip`. It runs PhiPack Profile in bounded
blocks, builds an uncorrected screening tree from the SKA alignment and reports
the root-to-tip permutation result. PhiPack is a detection screen in this mode;
it does not localise tracts or create a filtered alignment. Fast mode does not
produce a time-scaled tree or circulation/introduction interpretation.

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
