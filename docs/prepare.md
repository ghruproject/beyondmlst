# Independent preparation

`chronoclade prepare` imports frozen canonical profiles or validates and republishes
an existing prepared dataset. It writes `dataset.json`, typed categorical allele
matrices, inspection CSVs, complete imported lineage and metadata evidence,
`readiness.json`, `report.html` and a checksummed `prepare.json` stage result.
The complete output directory is published atomically; an existing output is never
replaced. Choose a new output directory for a revised input.

```console
chronoclade prepare profiles.json --catalogues catalogues.json --out prepared/
chronoclade prepare prepared/dataset.json --input-kind dataset --out imported/
chronoclade cgmlst prepared/dataset.json --out cgmlst/
```

Preparation is entirely offline in this first independent implementation. It
imports existing evidence without installing typing tools, querying reference
databases, resolving collection/accession identities or discovering public context.
Existing local assemblies explicitly linked by the input are copied into the
bundle under content hashes. The saved samples table uses relative assembly paths,
so moving the output does not break those references. Remote assembly URLs or
accession identities remain external references; no assembly download occurs.

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

## Readiness and remaining provider work

Dataset publication has a separate completion status from typing readiness. The
audit marks available profiles, absent loci, incomplete catalogue universes and
lineage assignment independently. Klebsiella readiness expects a resolved scoped
cgLIN assignment, including the authentic provider status `complete`; partial and
provisional assignments remain incomplete. E. coli expects resolved HierCC. Other species have lineage readiness
explicitly unassessed. A profile with missing lineage is incomplete. Missing
metadata is counted and retained; it does not exclude valid samples or by itself
make existing typing incomplete.

The audit says typing tools and reference database availability were not checked
because frozen imports do not need them. It never claims fresh typing was run or
that an authenticated reference database is ready. `prepare.json` lists relative
artifact paths, byte counts, hashes and exact sample IDs; input JSON and catalogue
snapshots are included as raw source evidence. Raw source values may retain their
original paths, while authoritative local assembly references are portable.

Collection membership resolution, exact accession identity/metadata retrieval,
input assembly acquisition for missing profiles, independent profile calling and
lineage assignment, and exact linked ENA/BioSample metadata enrichment remain to
be extracted behind query-only provider interfaces. The cumulative resolver mixes
those responsibilities with public context discovery and partition policy; the
independent prepare command deliberately does not wrap that runner.
