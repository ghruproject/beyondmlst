# Staged profile-first workflow

`chronoclade run` is now cumulative. Its default `fast` mode builds a profile
report; `full` adds selected context assemblies and the corrected genomic
analysis; `finish` runs the temporal tests and time tree when the clock gate
passes. Each completed stage keeps its own report and evidence files.

## Run a stage

The query dataset can be a metadata CSV, a Pathogenwatch collection, or a list
of accessions. For example:

```bash
pixi run chronoclade run queries.csv --mode fast --output results
pixi run chronoclade run --collection COLLECTION_UUID --mode fast --output results
pixi run chronoclade run --accessions accessions.txt --species "Klebsiella pneumoniae" --mode fast --output results
pixi run chronoclade run queries.csv --mode full --output results
pixi run chronoclade run queries.csv --mode finish --output results
```

The current options are shown by `pixi run chronoclade run --help`. Relevant
defaults are:

| Option | Default | Meaning |
| --- | ---: | --- |
| `--mode` | `fast` | Stop after the profile report (`fast`), corrected-tree stage (`full`), or dating (`finish`) |
| `--profile-limit` | `500` | Maximum public context records in profile analysis; query records are additional |
| `--context-size` | `50` | Maximum public context assemblies selected per lineage; focal query genomes are additional |
| `--nearest-per-query` | `3` | Nearest public profile candidates considered for each query during assembly selection |
| `--profile-bootstraps` | `10` | Requested locus-bootstrap replicates for descriptive group support |
| `--group-distance` | `0.02` | Exploratory complete-linkage cgMLST mismatch fraction |
| `--date-randomisations` | `100` | Tip-date permutations in `finish` |

`--collection` accepts a Pathogenwatch collection UUID or full collection URL.
Collection members are the query set. `--accessions` reads one identifier per
line, or a metadata CSV can supply accession fields. An accession list without
a frozen catalogue uses exact public identity matches; ambiguous or missing
matches stop with an error. `--catalogue` and `--public-typing` allow frozen
local reference inputs. `--query-typing` imports verified query assignments,
and `--typing-config` runs the configured native caller for query assemblies.

## What each mode produces

| Mode | Work and report |
| --- | --- |
| `fast` | Types available query assemblies when configured, resolves query and bounded public typing context, and writes `fast/profile_report.html` plus `fast/profile_analysis.json` and supporting results. It does not download context assemblies or build a recombination-corrected tree. The report shows input/profile coverage, exclusions, country figures, PCoA and neighbour-joining views where available, descriptive groups, exploratory allele-unit root-to-tip summaries, and nearest profile relatives. |
| `full` | Reuses the profile results to select up to 50 public context assemblies per lineage, with all focal query genomes retained, then runs the corrected genomic workflow. The stage landing page is `full.html`; lineage reports are saved as `report.full.html` beside the ordinary workflow outputs, with `supporting_results.full.zip`. This stage does not run date randomisations or create a time tree. |
| `finish` | Runs `full` and then the temporal tests. It writes `finish.html`, lineage `report.finish.html` files and `supporting_results.finish.zip`. A dated tree is produced only if the temporal evidence gate passes. |

The top-level `index.html` links to completed stages. `stages.json` records
their paths and the input fingerprint. A changed input invalidates saved finish
reports; the workflow checks the lineage-stage fingerprints before reusing
genomic results.

## Read the profile report

ChronoClade separates two descriptions of genetic groups. **Persistence across
years** reports the observed collection years represented by a group.
**Concentration in time and place** reports the observed year/location cells
and their share among records with both kinds of metadata. Both are descriptive
summaries of sampled records. Observations across years do not establish
uninterrupted persistence, and a concentrated cell does not establish an
outbreak or transmission link.

The report separates query and public-context record counts and profile
coverage. Records without a usable profile remain in input coverage and the
exclusion table. They cannot enter a comparable-profile cohort or its PCoA,
neighbour-joining tree, country figure or nearest-relative ranking. Country
figures describe profile-available members of a complete-comparability cohort
and combine query and context records; read them alongside the separate
query/context counts. The report also gives country and region summaries by
origin for all resolved input metadata, including records without a usable
profile, and a separate summary for the frozen public catalogue. These
denominators describe different pools and should not be conflated. Unknown
country and region values remain visible as explicit categories.

Each distance is a mismatch fraction over jointly called loci and reports the
number of shared called loci. For exports without a declared complete scheme
universe, comparisons use the observed locus union available to the compatible
records. The report warns when canonical scheme completeness is unknown; that
union must not be read as full scheme coverage. Incompatible scheme versions,
species, or locus sets are kept out of the same comparison cohort.

Genetic groups use deterministic complete linkage at `--group-distance` and
are candidate summaries for these data. Bootstrap support is reported only
when valid replicates and within-group pairs make it estimable. If no valid
replicates are available, support is unknown; that result is not evidence that
a group is unstable. A singleton has no within-group pair support to estimate.

PCoA and neighbour-joining trees summarize allele-profile distances. The
root-to-tip diagnostic uses an arbitrary neighbour-joining root and date
interval midpoints; it is exploratory and never passes the final clock gate.
The `finish` stage evaluates dates on the corrected genomic workflow and
creates a calendar tree only when the configured temporal screen supports it.

Location networks are possible maximum-parsimony state changes across alternate
sample roots. The root fraction records how often a possible change appears in
those tested rootings; it is not a probability or confidence interval. These
changes depend on the selected tree and the submitted locations, and do not
prove transmission, direction of spread or acquisition location.

## Current typing-data availability

The live request for the production Klebsiella cgMLST and LIN data is pending.
Without usable exported profiles, the workflow still retains query/context
metadata and reports missing-profile coverage, but profile distances, groups,
and their figures may be unavailable. `--public-typing` currently imports a
local frozen assignments file; it does not grant access to a protected live
database. Query assemblies can be typed only when a matching native caller and
reference database are configured with `--typing-config`. A user's
Pathogenwatch API key enables the collection analysis exports and authenticated
sequence downloads; without it, export provenance records the missing
credentials and the report shows the resulting coverage.

The corrected-tree stage also depends on retrievable assemblies for selected
context records. A record with a public typing assignment but no available
assembly can inform profile analysis while remaining outside the corrected
tree. The reports and selection audit keep profile-catalogue counts separate
from selected tree counts. A refreshed catalogue or different assembly
availability can change the selected comparisons and nearest relatives.
