# ChronoClade suite architecture review

Reviewed 10 October 2026 in the `codex/shared-suite-foundations` checkout.
This is a strict source-level architecture audit and implementation guide. The
five-stage design is in [tool-suite-design.md](tool-suite-design.md); this review
does not certify that those commands or their scientific workflows are complete.

## Verdict and scope

The reusable foundation is worth extracting now. The existing cumulative
fast/full/finish workflow is unsuitable as the permanent implementation of
prepare/cgmlst/esm2/tree/time. Adding command wrappers or more mode checks would
retain the central coupling, file-discovery semantics and shared invalidation
logic. Move ownership first, then give each stage its own validated inputs and
output namespace.

The review inspected the CLI, input/provider resolution, adaptive/refined context,
profile analysis/network/reporting, assembly/temporal orchestration and the ESM2
inference foundation. It also inspected current tests for selection, stage
snapshots, location reconstruction, country-tree display and widget interaction.
Line counts below are the baseline observed before concurrent foundation edits;
they are evidence of responsibility concentration, not a claim that this branch
created all the debt.

| Baseline module | Lines | Responsibilities observed |
| --- | ---: | --- |
| `report.py` | 1,414 | Theme, assets, downloads, context/scientific result sections, full/fast/summary pages |
| `profile_report.py` | 1,255 | Flexible dictionary decoding, tables, location views, profile/corrected reports and indexes |
| `profile_inputs.py` | 1,165 | Collection/accession parsing, identity joins, metadata overrides, exports, typing, context pool and assembly materialisation |
| `lineage.py` | 1,084 | Native execution/cache, reference/alignment/tree/recombination, temporal tests, dating, location inference and reports |
| `context_geography.py` | 996 | Identity deduplication, geographic summaries, plotting, focused views and output publication |
| `profile_analysis.py` | 892 | Profile compatibility/distances, groups/bootstrap, NJ/PCoA/regression, location outputs, table writes and display pruning |
| `staged_workflow.py` | 520 | Representative selection, cumulative stages, cross-stage invalidation, report snapshots and link rewriting |

Four existing files exceed 1,000 lines and two sit near it. No new suite feature
should be added into those files without first extracting the responsibility it
needs. Splitting a large function into neighbouring files while retaining every
mode flag is insufficient: the consuming boundary must become simpler.

## High-conviction findings

### 1. Cumulative orchestration prevents genuine independent stages

**Priority: P1 for suite completion.** `run_staged_workflow` resolves inputs and
runs profile analysis before every full/finish invocation. It derives one broad
fingerprint from profile/selection and temporal settings, validates selection
dictionaries inline, materialises assemblies, calls `run_workflow` in corrected
and full modes, then snapshots reports. `_run_lineage` in `lineage.py` branches
between corrected, fast and full behaviour and owns both tree construction and
dating. The CLI exposes the legacy runner rather than independent stage contracts.

The simplification is to publish `dataset.json`, partition/selection manifests,
`tree.json` plus its temporal assessment, and an independent time result.
Tree reads a saved selection; time reads a saved tree. Reports are written once
inside the owning stage's namespace. This removes orchestration-level scientific
mode checks rather than moving them into a new dispatcher.

`_snapshot_stage` is concrete deletion potential: recursive extension filters,
suffix/date exclusions, path relocation, HTML regex rewriting and copied bundles
exist to compensate for a shared mutable output directory. Keep existing snapshot
behaviour during migration, then remove it when independent stage output and
resume tests prove that it is no longer needed. Do not remove this protection
before replacement ownership works.

**Acceptance:** time-parameter changes leave tree and cgmlst artifacts unchanged;
two selections have separate outputs; failed time retains valid tree output;
tree/time startup does not invoke provider discovery or cgmlst analysis.

### 2. Location computation is reusable, but its current wrapper owns cgLIN policy

**Priority: P1 for reuse; extraction addressed by this foundation change.** Before
extraction, `profile_network.py` contained both generic fixed-tree
Sankoff reconstruction/exact pair extrema, graph metrics and layout, and
`build_profile_network`'s Klebsiella cgLIN-prefix subgroup routing. It imported
private `_compatible`/`_lineage` policy helpers from `context_refinement.py`.
`profile_country_tree.py` imported `_date_interval` from `profile_analysis.py` and
hard-coded NJ/cgMLST titles and branch units. Reusing those files directly in
tree/esm2 would import another stage's policy and mislabel scientific evidence.

Move fixed-tree reconstruction, centralities, colours, static network, interactive
widget and country-tree rendering into `location_network/`. Keep declared group
views outside the common engine. Renderers accept source tree kind, branch units,
root and analysed IDs. Date intervals belong in shared metadata handling.

Preserve the approved directed network, interaction, colour calculation, coherent
representative history and exact count ranges on unchanged fixtures. Do not
silently replace the scientific method during extraction. The older
`country_network.py` TreeTime marginal-state product has different uncertainty
semantics; its source parsing stays in a separate marginal adapter. Shared
presentation does not make marginal support and exact parsimony ranges equivalent.

ESM2 metadata composition and metadata-coloured embedding points can reuse these
visual components. Embedding coordinates or embedding NJ cannot be supplied as a
biological tree for ancestral location changes. An esm2 network needs an explicit
conventional cgMLST/SNP baseline, labelled with its source, or is unavailable.

**Acceptance:** unchanged fixed-tree edges/ranges/metrics/palette; no shared
location module imports profile_analysis/context-refinement stage policy;
corrected-SNP trees have SNP labels; embedding-space input is rejected at the
future typed biological-tree boundary. Filtering a saved graph cannot claim a
new reconstruction.

The current import graph meets the extraction boundary: common functions live
under `location_network/`; `profile_network.py` retains only the cgLIN view
adapter; shared country-tree rendering imports `metadata_dates.date_interval`
and receives scale/title/dating labels. `tree_basis` rejects named embedding
inputs and requires a Bio.Phylo tree, saving the declared basis in its result.
It is a caller declaration, not validation
of a frozen tree manifest or saved root/checksum provenance; those remain part
of the independent-stage contract.

### 3. Broad dictionary contracts hide identity, date and compatibility invariants

**Priority: P1 for correctness and maintainability.** `resolve_profile_inputs`
returns nested `queries`/`context`/`catalogue_rows`/`provenance` dictionaries.
`profile_report.py` repeatedly uses `_mapping`, `_items` and `_records` to accept
multiple result shapes and silently empty missing fields. Dataset identity,
profile/database scope, assembly paths and date eligibility are carried by
conventions rather than validated entry points.

Before extraction, date logic was duplicated across Pathogenwatch
`_date_bounds`/`normalize_dates`, profile analysis `_date_interval`, metadata
`_validate_date` and temporal serialization. They accept different formats:
the provider decoder can preserve date
intervals while the profile helper only parses year/month/day. That divergence
cannot be solved by a renderer choosing whichever string it recognises.
Identity/deduplication is likewise implemented across the provider catalogue,
cgLIN accession matching and geography aliases. These are not safely replaced by
one broad alias-set merge: ambiguous BioSample/assembly/run relationships need
typed evidence and explicit resolution status.

Make `datasets/` the owner of immutable IDs, typed tables, field provenance,
profile compatibility and versioned bundle validation. Preserve raw provider
values at the provider boundary. Standardise date intervals/statuses in one
metadata module, then adapt native backend formats explicitly. Downstream
analysis and rendering consume validated records rather than interpreting raw
provider variants. Stage-specific result payloads remain stage-specific.

**Acceptance:** portable round-trip preserves IDs, date precision/intervals,
missing calls, lineage scope and conflicts; ambiguous aliases cannot resolve by
first match; missing metadata does not remove valid genomes; malformed manifests
fail with a field-specific explanation rather than an empty report section.

The current foundation removes the provider dependency: `normalize_dates` and
the existing profile interval parser now live in `metadata_dates.py`, with
providers, datasets, profile analysis and fast-workflow callers importing their
canonical helper directly. It also moves scheme/version/frozen-export scope and
compatibility into neutral `typing_scopes.py`. `datasets/` imports no provider
or context-selection module. These moves establish ownership; reconciling every
accepted backend date format is still a separate task.

The review found and corrected an adapter compatibility blocker: records with
one scheme/version but conflicting nonempty database SHA256 values initially
entered a single matrix. `LocusCatalogue.database_sha256` and exact matching now
prevent that merge; known and unknown fingerprints also require separate explicit
catalogues. The original scope compatibility remains available to context
comparison through the neutral helper.

### 4. Selection is a shared policy and currently misses approved mandatory semantics

**Priority: P1 before tree-contract acceptance.** Before extraction,
`select_assembly_context` lived inside the orchestration module, while fast
workflow imported it back from that module. Both now use `selections.py` as the
canonical owner; exact pool/decision validation is also owned there. Its distance
reader still reconstructs a full in-memory map from pairwise CSV.
It slices each query's candidates to `nearest_per_query`, then applies a capped
nearest-neighbour budget. Pins above the context budget raise rather than being
retained with an explicit overrun. These behaviours do not satisfy the design's
preserved ties and mandatory-neighbour/pin rules.

The move into `selections.py` fixes ownership without changing selection results.
Later replace the capped queue with a mandatory-set union:
inputs, pins and declared required nearest groups/ties first; quota/diversity
selection fills the remaining budget. A mandatory overrun is recorded, never
silently truncated. Pass a distance-evidence reader/access interface rather than
making selection depend on profile_report paths or a CSV-only representation.
Do not turn selection into a shared distance engine.

**Acceptance:** same existing selections after extraction; separate new fixtures
for exact ties, overrun, missing dates, reproducible alternatives and sources from
cgmlst/esm2. Every selection has source IDs/hashes and per-ID reasons. Tree checks
that emitted manifest and uses exactly its IDs.

### 5. Reporting crosses private boundaries and owns some scientific assessment

**Priority: P2.** `profile_report.py` imports private geography/context/network/
evidence/recombination helpers from `report.py` alongside its styles. This means
editing the full report can break profile reports and prevents a clear report
owner. `temporal_report.py` also contains `assess_temporal_signal`, an assessment
used by scientific orchestration, alongside fonts and plot writers.

The theme/fonts now live in `report_components/styles.py`, imported directly by
full/profile/stage report writers. Follow with public components for figures, tables, downloads and
named result sections. Move temporal assessment into the scientific temporal
layer; report writers display its saved result. Stage report functions retain
their result-specific content rather than one generic renderer with mode checks.
Avoid a universal recursive dictionary normaliser: explicit adapters at the
old/new boundary are easier to delete and inspect.

**Acceptance:** profile/corrected/full reports retain approved layout and labels;
new shared report components never call scientific engines or provider APIs;
full and profile writers no longer import each other's private helpers;
report rebuild works from saved evidence without rerunning analysis.

### 6. Execution and hashing have several competing owners

**Priority: P2 before independent resume is claimed.** `_run_command` and direct
file fingerprinting live in `lineage.py`; temporal workers execute commands
separately; stable object hashing is imported from the Pathogenwatch provider;
file SHA256 exists independently in temporal, ESM2 proteins and neighbourhood.
JSON publication is atomic in some ESM2 paths while legacy stage/report manifests
are written directly. The cumulative fingerprint mixes profile and temporal
settings, so changing dating options can invalidate reader-facing tree products.

Create one execution/artifact layer when migrating independent stages. It owns
content/file hashing, validated relative paths, tool/version-aware job identity,
atomic publication and a resource plan. The provider client must not be the
canonical hashing utility. Reuse one ordinary job description under local and
SLURM executors; do not implement two scientific pipelines.

Preserve the useful ESM2 atomic publication and provenance-specific vector-cache
behaviour. Backend/device fallback must not mix caches. Do not force the model
cache into the same key schema as a native alignment job simply because both
contain hashes.

**Acceptance:** changed direct input or tool identity invalidates its owner;
unchanged upstream bundles remain valid; interrupted publication preserves the
previous complete manifest; unavailable explicit device/site resource fails
clearly; independent work remains bounded by the shared CPU allocation.

### 7. Scale and optional-dependency claims need separate completion gates

**Priority: P1 for large-block claims, P2 for environment separation.**
`profile_analysis.py` enforces `MAX_RECORDS = 1500`, builds dense pairwise matrices,
uses full eigendecomposition and Python NJ, and writes every pair to CSV. Its
greedy complete-comparability cohorts are reported limitations, not proof of
biological partitioning. Removing the guard would expose resource failure without
making the engine a tens-of-thousands implementation.

The ESM2 foundation is well scoped: `esm2_cli.py` lazily imports protein inference,
the engine retains sequence mapping/model/device/cache provenance, and describes
itself as protein-level only. It is not the prepared-dataset genome embedding,
selection or report route. Keep that distinction. `pixi.toml` still installs all
native tools in its default environment; `cli.py` imports workflow/lineage/report
eagerly. A base CLI that avoids PyTorch is useful but does not complete all stage
environment boundaries.

Replace dense distance/tree/ordination/storage backends behind validated small-
input interfaces and benchmark complete outputs. Split tool environments and lazy
stage imports independently. Optional model inference can progress without being
used to imply biological validation or suite integration.

**Acceptance:** unchanged small-input distance evidence; full IDs retained at
benchmark sizes; measured memory/runtime for network reconstruction too; CLI and
cgmlst work without PyTorch/weights/tree/date executables; model inference and
actual genome-level neighbour/selection performance are tested separately.

## Foundation change and follow-through

The shared branch starts by extracting real reused code into `datasets/`,
`location_network/`, `metadata_dates.py`, `selections.py` and
`report_components/styles.py`, with neutral `typing_scopes.py` for namespace
identity. These are reviewable ownership improvements.
The old `country_network.py`, `profile_country_tree.py` and
`profile_network_widget.py` implementations are removed, with callers/tests
migrated to their shared owners; no compatibility forwarding modules remain.
The corresponding network/date/style/selection migrations should preserve
behaviour, and any new dataset contract should clearly state its schema scope.
Passing their focused tests does not certify all findings above as resolved.

An independent AST comparison against the pre-extraction source found unchanged
Sankoff reconstruction, root-pair sensitivity, centralities, nearest-observation
aggregation, static layout/drawing, marginal network calculation/reporting,
representative selection/distance reading and profile date parsing, apart from
renamed function identifiers and formatting. The shared renderer now accepts
method-specific text; profile call sites explicitly pass the original labels.
This is evidence of preserved calculation, not a substitute for the branch's
focused tests and rendered-report checks.

The final independent foundation review exercised dataset, shared-location and
shared-component tests: 51 passed. It rechecked provider/stage import boundaries,
the conflicting/unknown database-fingerprint regression fixtures and declared
tree-basis persistence. No remaining blocker was found within this extraction
and prepared-bundle scope. This approval does not cover independent prepare/
cgmlst/tree/time commands, the revised mandatory-selection policy, scale
backends or genome-level ESM2 analysis.

Remaining work should follow the parallel packages and integration gates in the
design document. The first complete vertical path is prepared bundle → explicit
cgmlst partition → selection manifest → independent corrected tree/assessment →
independent dating result. Then exercise alternatives, resume/invalidation,
large-block execution and the experimental ESM2 route separately.

Do not combine a scientific algorithm change, a whole-report restyle and an
ownership extraction in one unreviewable diff. Maintain fixed-tree numerical
fixtures and rendered reports as comparison evidence. The approval bar is a
smaller dependency graph and fewer coupled responsibilities, alongside preserved
scientific meaning and visible output, not merely passing legacy mode tests.
