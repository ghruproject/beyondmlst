# Frozen public Pathogenwatch metadata

Pathogenwatch context starts with a complete public same-ST metadata snapshot,
before QC filtering, focal exclusions, dated-cohort eligibility or downsampling.
A catalogue describes available sequenced records; it does not estimate country
prevalence, incidence, migration direction or transmission.

## Contract and public scope

The [Pathogenwatch API documentation](https://pathogen.watch/docs/api) embeds its
OpenAPI specification in the page's `__NEXT_DATA__` JSON. Verified on 9 October
2026, the supported-organism endpoint is `GET /api/organisms/supported`. It lists
MLST and organism-specific analysis capabilities separately. Genome availability
or cgMLST support does not imply cgLIN support; `LIN Codes` is a separate
Klebsiella pneumoniae capability.

The provider also supports *E. coli* (organism `562`), whose advertised analyses
include `HierCC`. Public genome details do not expose allele profiles or HierCC
codes directly. Use a frozen `--public-typing` import or native typing of the
bounded downloaded pool, rather than assuming an undocumented HierCC download
job. `--scheme ecoli` queries the primary `mlst` field. See
[new assembly typing and HierCC](query-typing.md) for refinement.

`POST /api/search/genomes` takes an explicit `organismId`, the ST array in the
explicit `mlst` or `mlst2` field, and `qc: [true, false]`. Pagination uses query
parameters `limit`, `sort=id` and `after=meta.endCursor`, with the total in
`meta.count`. The documentation does not specify a search page-size ceiling or
a snapshot-isolation guarantee. The client defaults to 100 records per page and
rejects count changes, missing/repeated cursors, repeated IDs, missing public
access labels and totals that do not reconcile. Empty results are valid.

Each search row's `uuid`, rather than its numeric `id`, joins
`GET /api/genomes/details?id=UUID`. Details must match UUID, organism and ST.
All records must be labelled `projectAccess=PUBLIC` in search; explicitly
non-public details are also rejected. These reads deliberately omit API keys
and session cookies, including when a download credential is configured. This
prevents an authenticated search from adding private genomes. Ordinary offline
CI uses injected transports and never requires credentials or the live service.

The live ST147 query on this date returned a reported count of 7,807. This is a
search observation, not a statement about a completed pilot catalogue or its
metadata completeness. Measure completeness from a successfully frozen snapshot.
Public details observed both null top-level location and nested `Country=IN`,
and metadata represented as an array of name/value pairs. The API schema is less
specific than the response; offline fixtures assert this observed contract.

`X-API-Key` is the documented authentication header for protected operations.
`load_api_key()` reads `PATHOGENWATCH_API_KEY`, or
`~/.config/chronoclade/pathogenwatch.json` with an `api_key` string. The file must
have user-only permissions (`0600`). No credentials enter catalogue provenance,
raw responses, error messages or query logs. Download/export code must separately
request authentication; public discovery never does.

## Python interface and snapshot schema

```python
from chronoclade.pathogenwatch import PathogenwatchClient, load_catalogue

client = PathogenwatchClient(workers=6)
snapshot = client.freeze_catalogue(
    "st147_catalogue.json", organism_id="573", st="147", mlst_scheme="mlst"
)
snapshot = load_catalogue("st147_catalogue.json")  # verified offline replay
rows = snapshot["rows"]
```

`search()` returns the exact query, raw pages, all search records, expected count
and explicit completeness. `freeze_catalogue()` fetches all detailed records
with bounded concurrent requests, retries transient HTTP/transport failures and
writes the completed JSON atomically. Non-transient failures, access denial and
contract mismatches raise `PathogenwatchError`; a partial acquisition is never
written as a complete catalogue. Failed freezes must be retried, since upstream
search does not promise a transactionally stable snapshot.

The envelope has `schema_version=1`, `rows`, `raw_search_pages`, `raw_details`
keyed by UUID, `raw_supported_organisms`, `deduplication_audit` and `provenance`.
Provenance retains exact query/hash, UTC retrieval time, raw count, sample-unit
count, page/detail counts, capability response, measured metadata completeness,
elapsed time and SHA-256 hashes of raw responses and normalised rows.
Unavailable service/database versions are explicitly `unavailable`.
`snapshot_sha256` covers the entire envelope excluding that field. The loader
checks these hashes/counts, public completeness and deterministic normalisation
of all raw records before returning the envelope. Snapshots are metadata, not
assembly caches, and contain no API credentials.

Each normalised row keeps `source=pathogenwatch`, `source_genome_id` (UUID),
`numeric_source_id`, source project/checksum, `sample_id`, species, organism,
MLST field and ST separately. `biosample`, `biosample_accessions`,
`run_accessions`, `assembly_accessions` and `study_accessions` retain recognised
accessions from explicit detail fields, names and metadata. Missing linkage stays
missing. Names and same ST do not prove biological identity. cgLIN fields are
added by the separate annotation module; no MLST value is borrowed as a cgLIN.

Country precedence is named nested metadata (`Country`, `country`, country of
origin, `geo_loc_name`), detailed location, then search location. `country_raw`
retains every input with its field name; `country_provenance` records the chosen
field and `country_conflict` flags disagreeing normalised countries. ISO codes
use frozen English CLDR territory labels; common aliases are merged, ENA place
suffixes after a colon are removed, and unmatched names are retained rather than
invented. Missing values become `Unknown`. Country is not country of acquisition.

`collection_date` preserves original precision and `date_start`, `date_end`,
`date_precision` retain year/month/day/interval bounds. Bounds are not asserted
exact collection days. Metadata dates take precedence over service date bounds;
when metadata is absent, whole-year/month bounds are recognised as such.
Invalid/reversed intervals are ineligible; undated rows stay in the catalogue.
`dated_cohort_eligible` only means a valid available date, not a passed temporal
signal gate. Host and isolation source are optional and never implicit exclusions.

`source_qc`/`qc_pass` preserve the provider's boolean QC status or `None`, and
`source_length`, `source_n50`, `source_contigs` preserve available metrics. These
are Pathogenwatch fields; they are not completeness or contamination scores.

## Sample counting and focal aliases

`deduplicate_catalogue(rows, focal_aliases=())` returns normalised sample units
and an audit. It forms deterministic connected components using BioSamples,
read-run accessions and version-independent assembly accessions, then prefers a
QC-pass representative and the lexicographically first source UUID. Shared
study IDs never link samples. Source UUIDs link repeated source records. It
retains `sample_unit_id`, all `source_genome_ids`, `raw_genome_count`,
`identity_resolved` and `identity_conflict`; multiple BioSamples in one component
are flagged as identity conflicts and are not called unique biological samples.
Unlinked records are separate source units. These units are not unique patients
or infections.

Focal overlap uses only supplied source UUIDs or strong sample/run/assembly
aliases, excluding the whole component with a recorded match. Each component
records its representative, source IDs, duplicate/overlap reason and metadata
conflicts. Pass focal accessions explicitly; do not exclude by human-readable
sample names or studies. Apply the intended QC rule before deduplicating for the
eligible geography denominator, while retaining the original envelope's raw
count and separate audit of losses. Dated selection is a subsequent filter and
must not remove undated country-summary records.

## Testing

`pytest tests/test_pathogenwatch.py` exercises multiple search pages, anonymous
public scope, detail joins, duplicate/focal aliases, country conflicts, date
precision, QC failures, missing metadata, transient failures/access denial,
incomplete pages, deterministic replay, credential permissions and corruption.
An opt-in live smoke should call `request_json` for a single search page with
`limit=2`, then inspect a returned public detail; `search()` traverses all pages.
Freeze the whole catalogue when full completeness metrics are required. Live results and large raw snapshots should remain outside
committed fixtures.
