# Country composition within cgLIN groups

`chronoclade.context_geography.generate_context_geography` regenerates figures and
source tables offline from the frozen, annotated normalised catalogue. It does not
require assemblies, phylogeny, dates or Pathogenwatch access. The workflow can embed
its returned `report_html` fragment and include `outputs` in supporting results.

```python
from chronoclade.context_geography import generate_context_geography

geography = generate_context_geography(
    frozen_annotated_rows,
    "results/context-geography",
    selected_source_ids=selected_pathogenwatch_ids,
    focal_rows=validated_annotated_focal_rows,
    scope={
        "description": "Public Pathogenwatch K. pneumoniae species-complex ST147 catalogue",
        "snapshot": "retrieval timestamp and frozen snapshot hash",
        "filters": "Explicit Pathogenwatch QC pass; deduplicated; undated records retained",
    },
)
```

Use the **full explicitly QC-passing same-ST catalogue before dated-cohort exclusions,
country/year balancing and assembly down-selection**. Restrict the catalogue to the
intended species/ST and public source before calling the function, and describe those
restrictions in `scope`. Unknown QC is excluded from this explicit QC-pass denominator.
The separate selected-context denominator follows the requested source IDs, including
all aliases of a deduplicated record. Focal rows have their own denominator; strong
accession overlap with public records is audited and never added to public proportions.
Focal dates or QC do not silently remove survey records from their separate panel.

The counting unit is a biological sample where BioSample/accession evidence supports
it, otherwise a record unit of unresolved identity. It is not necessarily a patient or
an infection. Connected BioSample, assembly and run aliases identify duplicates;
`sample_unit_id`, `source_genome_ids`, `raw_genome_count` and `identity_resolved` from
the provider deduplication are supported. The JSON audit retains raw counts,
representative IDs, raw/normalised countries, identity uncertainty, conflicting sample
identities, country conflicts and focal overlaps. Conflicting country values within a
sample unit become Unknown rather than borrowing one record's country. Country
normalisation and raw-value/provenance retention occur in the provider normaliser.

cgLIN annotations use `cglin_group_5`, `_6` and `_7` and their corresponding
`cglin_status_<depth>` fields. They must be validated full-prefix keys scoped by
scheme and version, produced by the cgLIN annotation importer. A partial code resolves
only its actual prefix depths. Conflicting duplicate assignments become unresolved.
Unresolved/missing/malformed/provisional-at-unresolved-depth statuses appear as
**assignment coverage**, outside biological lineage groups. Resolved provisional
prefixes remain explicitly marked provisional in tables. Scheme/version and full
prefix appear in lineage labels; equal final components never merge different groups.
Depths 5, 6 and 7 are exploration levels, **not pairwise SNP cutoffs**. Group-size,
singleton and assignment-coverage summaries are returned for evidence-based default
choices; this module does not pick a universal preferred depth.

`country_composition.csv` and `.tsv` retain every named country, every singleton and
small group. Each row includes source, snapshot, scope, filters, scheme/version, prefix
depth/key, assignment category/status, cohort, country, sample-unit count, raw record
count, total denominator, known-country N, Unknown N and percentage. Percentages are
`100*n/N`, including Unknown, and sum to 100% within each group subject to floating
point rounding. The raw-record counts reconcile with the eligible raw catalogue while
the deduplicated counts reconcile with group plus assignment-coverage totals.

Horizontal count and percentage plots export publication-sized SVG and 180 dpi PNG.
Groups paginate at 24 per figure. The visual displays the top 11 named countries by
full public-catalogue sample count, plus Unknown; remaining countries form Other.
This fixed public-country rule applies across depths. Bounded selected/focal cohorts
show all named countries when at most 12 are present; their captions expose the rule.
The full table remains available.
Country colours are deterministic across runs and panels; Unknown and Other have
distinct reserved colours. Figures label sample-unit denominators, known/Unknown N,
singletons and small groups. Empty cohorts or unsupported/missing annotation show a
clear message; no lineage is invented from unavailable codes. Real unresolved records
can still have a separate assignment-coverage country figure.

These figures describe **the geographical composition of available sequenced records**,
affected by surveillance/submission coverage and selection. They do not estimate
country prevalence, incidence, migration direction or transmission. Country metadata
is not country of acquisition. The same-ST catalogue is not automatically the entire
global cgLIN group across all STs. Retain these captions when using exports elsewhere.

Offline tests use hand-counted multiple-country groups, prefixes sharing their final
component, Unknown country, partial codes, duplicate BioSamples/assemblies, focal
public overlap, singleton groups and deliberately imbalanced selection. The release
pilot should additionally inspect its real SVG/PNG and HTML at readable size, verify
counts against frozen raw exports and archive the source tables and snapshot hashes.
