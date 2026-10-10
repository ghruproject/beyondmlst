# Independent preparation

`chronoclade prepare` resolves live query collections, accessions or assembly
metadata, imports frozen canonical profiles, or validates and republishes an
existing prepared dataset. It writes `dataset.json`, typed categorical allele
matrices, inspection CSVs, complete imported lineage and metadata evidence,
`readiness.json`, `report.html` and a checksummed `prepare.json` stage result.
The complete output directory is published atomically; an existing output is never
replaced. Choose a new output directory for a revised input.

```console
chronoclade prepare profiles.json --catalogues catalogues.json --out prepared/
chronoclade prepare prepared/dataset.json --input-kind dataset --out imported/
chronoclade cgmlst prepared/dataset.json --out cgmlst/
```

Frozen profile and dataset imports are offline. Live preparation uses the same
portable bundle contract, with exact query identity, provider response snapshots,
metadata precedence/conflicts and per-sample retrieval/readiness evidence. It does
not discover public context or download context assemblies. Linked local input
assemblies are copied under content hashes and remain portable.

## Live query preparation

```console
chronoclade prepare https://pathogen.watch/collections/jX5cwsoyJ1KDquMqssUzAD \
  --input-kind collection --out prepared-collection/
chronoclade prepare accessions.csv --input-kind accessions \
  --species "Klebsiella pneumoniae" --metadata overrides.csv --out prepared-accessions/
chronoclade prepare assemblies.csv --input-kind assemblies \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json --out prepared-assemblies/
```

Collection identifiers accept a short UUID or a full HTTPS Pathogenwatch collection
URL, including its optional slug. Membership pages must reconcile exactly with the
advertised count. An accession input is a unique one-ID-per-line list or CSV with
`accession` or `source_genome_id`; CSV metadata is applied as an explicit override.
Public accession search requires species, a complete bounded response and exactly
one exact accession match. Ambiguous identities fail without choosing a neighbour.
An assembly input is a CSV with `sample_id`, `assembly` and `species` (or
`--species`). Assembly paths are relative to the CSV. Additional `--metadata`
overrides join by exact ID and take precedence over provider values; differences
remain in the conflicts table.

Existing server-advertised cgMLST and cgLIN analysis exports are retrieved for the
resolved input IDs. The configured Pathogenwatch key is read from the environment
or the protected local configuration; it is never accepted as a command-line
argument or written into snapshots. Public searches remain credential-free.
`--query-typing`, `--public-typing` and `--cglin-export` accept validated frozen
assignments. Query typing must match exact sample identity, species and any known
input assembly SHA256. Missing profiles and lineages remain separate checks.

`--typing-config` enables the pinned native upstream caller and assigner. A
compatible existing profile missing its lineage goes directly to the assigner,
without assembly acquisition. Other incomplete inputs acquire **input** assemblies
and run native typing. A missing/unready configuration publishes an incomplete
bundle with an explicit reference-database reason, rather than inventing calls.
Typing-tool installation and database readiness are separate audit fields. Existing
frozen typing does not require a Pasteur key. Fresh Klebsiella lineage assignment
requires its configured Pasteur reference database; E. coli assignment requires
its scheme-specific HierCC reference database.

Live metadata enrichment uses exact linked ENA run/BioSample accessions. It fills
missing country, collection date, host and isolation source only when returned
records agree. Existing values retain precedence, source raw values and conflicts
remain recorded, and unrelated identities are rejected. It never derives a region
or date precision from an upload date. `--no-enrich-metadata` disables enrichment.
Missing metadata and request failures preserve the sample.

Provider cgMLST exports list observed loci only. Their automatically generated
catalogue is explicitly incomplete; no full-universe call fraction is claimed.
Supply `--catalogues` with a matching full scheme catalogue to assess completeness.
Native callers provide the ordered full locus universe. Source responses and exports
are copied into the checksummed bundle after credential checks.

## Frozen input contract

Supply a JSON list of canonical records, or an object with `records` plus optional
`lineages`, `crosswalk`, `provenance`, `conflicts`, `exclusions`, `retrieval` and
`parameters`. Evidence tables are lists; parameters is an object. Explicit evidence
tables override the frozen-record adapter's automatic evidence. Every evidence row
must refer to an existing `sample_id`. Unknown wrapper fields are rejected. The
legacy resolver's `queries`/`context` output requires explicit conversion to this
contract; the new command does not run the legacy resolver.

```json
[
  {
    "sample_id": "input-001",
    "role": "input",
    "species": "Klebsiella pneumoniae",
    "mlst_st": "39",
    "collection_date": "2024-02",
    "country": "Greece",
    "cgmlst_scheme": "example:cgmlst",
    "cgmlst_scheme_version": "v1",
    "cgmlst_status": "resolved",
    "cgmlst_profile": {"locus-a": "001", "locus-b": "7"},
    "cglin_scheme": "example:cglin",
    "cglin_scheme_version": "lin-v1",
    "cglin_status": "resolved",
    "cglin_raw": "1.2.3.4.5.6.7"
  }
]
```

Records require safe unique `sample_id` values. `role` is `input` or `context`;
the existing `origin` values `local`, `query`, `focal`, `input`, `context` are
also accepted. A record with neither role nor origin defaults to `input`.
Explicitly supplied context records are preserved without discovery or selection.
`--species` fills missing species values and rejects conflicting values. It cannot
override an existing prepared bundle's species. Date precision, null metadata,
ambiguous identities and complete lineage assignment strings/maps remain saved.
Profile allele values are categorical; missing calls never become allele matches.

An explicit catalogue JSON list is required for records, even when profiles are
unavailable. An object containing only `catalogues` is also accepted.

```json
[
  {
    "scheme_id": "example:cgmlst",
    "scheme_version": "v1",
    "loci": ["locus-a", "locus-b", "locus-never-called"],
    "source": "frozen-schema.json",
    "complete": true,
    "database_version": null,
    "database_sha256": null
  }
]
```

`complete` must be explicitly `true` or `false`. Do not claim completeness from
a union of observed calls. Every declared locus is retained, including loci absent
from every sample. Catalogue scheme/version and any database SHA256 must match
the records; incompatible called profiles fail rather than being discarded.
Missing profiles remain as samples, including all-null matrix rows when their
scope matches a supplied catalogue. Credentials are rejected before source
snapshots or outputs are published.

## Readiness and public profile context

Dataset publication has a separate completion status from typing readiness. The
audit marks available profiles, absent loci, incomplete catalogue universes and
lineage assignment independently. Klebsiella readiness expects a resolved scoped
cgLIN assignment, including the authentic provider status `complete`; partial and
provisional assignments remain incomplete. E. coli expects resolved HierCC. Other species have lineage readiness
explicitly unassessed. A profile with missing lineage is incomplete. Missing
metadata is counted and retained; it does not exclude valid samples or by itself
make existing typing incomplete.

Frozen-import audits say tools/databases were not checked because those imports
do not need them. Live audits distinguish configured, unavailable and ready
references, and whether native typing or assignment actually ran. Dataset
publication completion does not imply typing completeness.

The cgMLST provider hook `discover_profile_context` consumes a prepared bundle and
publishes a new prepared bundle for profile analysis. It searches the complete
accessible public same-ST pool, records counts at LIN5/6/7, and applies the explicit
chosen LIN level (default 5). E. coli requires an explicit HierCC level. Lineage
scheme/version/database scope and profile scope must be compatible. The hook
requests every matching profile, without a candidate limit or context assembly
acquisition. It deduplicates source identity, preserving different samples with
identical profiles; exact inputs are excluded. Unavailable or incompatible profiles
remain audited exclusions.

`readiness.json` and dataset parameters record discovered, requested, retrieved
and usable counts separately, pagination evidence, export status, exclusions,
source dataset SHA256 and the explicit policy. Input metadata, profile calls and
lineage/crosswalk/retrieval evidence are preserved. Accessible public matches are
reported; global completeness is never claimed. Raw provider snapshots are
checksummed under `sources/provider/`.

On 10 October 2026, live query-only preparation retrieved all 23 supplied
*S. aureus* collection members and all 23 existing cgMLST profiles without assembly
acquisition or context discovery. The observed locus catalogue remained incomplete
and lineage readiness was explicitly unassessed for this species. This validates
live membership/export acquisition, not Klebsiella native reference-database
readiness. A live public Klebsiella ST147 search reconciled 7,807 accessible
records across 79 pages; this verifies search pagination, not complete lineage or
profile export retrieval. Local upstream MLST 8.0.0, plincer 7.0.0 and hclink 4.0.1 installations
were discovered; their reference databases were unavailable, independently of
those installed tools. A fixture with 120 matching public profiles verifies the
unbounded explicit-level context path and exclusion of an unrelated prefix.

Upstream tool contracts are documented by [Pathogenwatch MLST](https://github.com/pathogenwatch-oss/mlst),
[plincer](https://github.com/pathogenwatch-oss/plincer) and
[hclink](https://github.com/pathogenwatch-oss/hclink). Pathogenwatch service endpoint
response shapes are frozen and validated by the adapter and its live checks;
provider response compatibility remains explicit rather than assuming a database
release from an analysis job name.
