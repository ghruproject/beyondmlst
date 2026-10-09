# Running an analysis

## Choose the stopping stage

`chronoclade run` defaults to `fast`, the profile-first analysis. Use `full` to
add selected context assemblies and the corrected genomic analysis. Use
`finish` to run those same steps and then assess temporal signal and create a
time-scaled tree only when the temporal gate passes.

```bash
# Profile analysis only; no context assemblies are downloaded.
pixi run chronoclade run metadata.csv --mode fast --output results

# Corrected genomic analysis; stops before date tests.
pixi run chronoclade run metadata.csv --mode full --output results

# Corrected analysis followed by temporal tests.
pixi run chronoclade run metadata.csv --mode finish --output results
```

Use `pixi run chronoclade run --help` for every current option. The query source
can be a metadata CSV, `--collection COLLECTION_UUID`, or
`--accessions accession-list.txt`. Accession lookup may require `--species` or
a frozen `--catalogue`. `--public-typing` and `--query-typing` load local
assignments; `--typing-config` runs the configured native caller for query
assemblies.

`--dry-run` writes `plan.json` and prints the planned stages. It does not fetch
or validate live typing data, download context assemblies, or run the analysis.

## Stage sizes and selection

The defaults are 500 public records for profile analysis, up to 50 context
assemblies per lineage in `full` and `finish`, three nearest profile relatives
per query considered during context selection, ten requested locus-bootstrap
replicates and a complete-linkage mismatch fraction of 0.02. Query records are
additional to both public-context limits. Override these with
`--profile-limit`, `--context-size`, `--nearest-per-query`,
`--profile-bootstraps` and `--group-distance`.

Context selection is limited to the recorded public profile catalogue. It
considers nearest compatible profiles, then representatives across genetic
group, year and region; remaining slots favour profile diversity where
distances are available, with metadata-only background records used when no
comparable profile exists. The selection audit records each reason. A selected
record still needs a retrievable assembly to enter the corrected tree.

## Temporal analysis

`--date-randomisations` defaults to 100, and applies in `finish`. One hundred
permutations give a minimum corrected empirical p-value of 1 / 101 = 0.0099.
The default `--date-randomisation-method root-to-tip` permutes dates on the
corrected tree; `full-tree` refits the complete TreeTime model and reroots on
each permutation. The date gate also requires a positive estimated rate and
completion of every requested permutation. R² alone does not pass the gate.

When fewer than three usable distinct collection dates are available, the
corrected report remains available and `finish` records that temporal signal
was not assessed. It does not produce a failed randomisation result. For
supported datasets, the time tree estimates ancestral dates under the selected
clock model; it is not a transmission tree.

## Outputs

Open `results/index.html` for links to completed stage reports. Fast mode writes
`results/fast/profile_report.html`, profile tables, figures and a supporting
results archive. Full mode adds `results/full.html` and one
`report.full.html` plus `supporting_results.full.zip` per analysed lineage.
Finish mode adds `results/finish.html` and `report.finish.html` plus
`supporting_results.finish.zip` for each lineage. Earlier stage reports remain
available after later stages complete.

The output contains `stages.json` with stage paths and a fingerprint of the
resolved query/context inputs. Reuse is based on those inputs and on each
genomic stage's direct inputs. `--force` requests rerunning completed genomic
stages.

## Existing context manifests

`prepare-context` remains available for users who need to inspect or freeze a
context pool before the staged run. It writes `combined_metadata.csv` and
`context_manifest.tsv`; the main command accepts that manifest with
`--context-manifest`.

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme klebsiella --st 147 --output context/ST147 \
  --candidate-pool 500 --max-context 50 --dry-run

pixi run chronoclade run context/ST147/combined_metadata.csv \
  --context-manifest context/ST147/context_manifest.tsv \
  --mode finish --output results
```

The profile-first and context-manifest routes both keep the scope of the public
search visible. A nearest relative is nearest only within the candidates that
were available to that run.
