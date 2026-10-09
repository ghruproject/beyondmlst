# Comparable cgLIN assignments

ChronoClade supports the **Klebsiella pneumoniae species complex
scgMLST629_S** nomenclature mirrored by Pathogenwatch from BIGSdb-Pasteur.
An available genome or an MLST ST does not establish cgLIN capability.
Other organisms and schemes receive an explicit unsupported status.

## Meaning and indexing

Depth is the number of components, counted from one. Depths 5, 6 and 7
are exploratory nested prefix groups, not pairwise SNP thresholds or
guaranteed transmission clusters. The full prefix defines a group: e.g.
`0,0,197,0,4` differs from `1,1,198,0,4`, despite the shared last component.
Machine keys encode `[scheme, scheme_version, [full prefix components]]` in
canonical JSON. Keep these keys separate from source codes and existing
ChronoClade species/ST lineage labels.

The scheme uses ten components. Complete numeric codes resolve all ten;
partial codes resolve only their observed consecutive numeric prefix.
Unknown suffixes (`-`, `?`, `*`, or provisional `*identifier`) never become
numeric assignments. Numeric values after an unknown component, negative
values, arbitrary text or more than ten components are malformed. Comma
separation matches the verified Pathogenwatch export; underscore separation
matches BIGSdb displays. Both normalise to identical numeric prefixes.

BIGSdb-Pasteur describes provisional cgSTs marked with `*`, and incomplete
codes derived from the shared prefix with an existing reference. ChronoClade
preserves that provisional flag and the code's independent complete/partial
status. It never expands a partial code using the closest profile field.
The verified CSV also represents provisional cgSTs as 40-character
hexadecimal profile hashes without the UI's asterisk; these remain raw and
are explicitly flagged, including when no LIN code is available.
See [BIGSdb-Pasteur cgMLST/LIN documentation](https://bigsdb.web.pasteur.fr/klebsiella/cgmlst-lincodes/).

## Verified export contract

The [Pathogenwatch API documentation](https://pathogen.watch/docs/api),
including its embedded OpenAPI specification, was inspected on 9 October
2026. The generic analysis download uses:

```text
POST /api/downloads/klebsiella-lincodes?job=lincodes-3390273-2
Content-Type: application/json
X-API-Key: <protected credential>
{"ids":"174724,174727"}
```

`job` is an existing **analysis job name**, not a new export job ID. The
live checked two-record public ST147 request returned HTTP 302 pointing to
an S3 CSV. The response was gzip encoded; the client decodes it before
freezing exact CSV bytes. Request IDs are numeric database IDs from the
frozen catalogue's `numeric_source_id`; export `Genome ID` values are UUIDs
and join to `source_genome_id`. A request using UUIDs produced a server
error, so the two identifier types must stay distinct.

The checked CSV headers were `Genome ID`, `Genome Name`, `cgST`,
`Closest cgST`, `LIN code`, `Sublineage`, `Clonal Group`, `Identity`,
`Identical`, `Compared Loci Count`, `Closest profile(s)`. Both public records
`tAYgKfJMnCBK81yvS4eiYN` / `SAMN20121908` and
`qAzBoQ6TzMVYhFLR5z6aq9` / `SAMN20121911` returned cgST `19550` and
`0,0,197,0,4,0,93,0,0,0`. Their depth 5–7 prefixes are
`0,0,197,0,4`, `0,0,197,0,4,0`, and `0,0,197,0,4,0,93`.

The public documentation does not specify polling, export URL expiry or
the server's maximum batch size. No polling protocol or expiry is assumed.
The implementation uses conservative sequential batches of at most 100
records; this is a local limit, not a claimed upstream limit. The checked
analysis job version is retained separately from the nomenclature version;
it is **not** evidence of a BIGSdb database release. Where the source gives
no nomenclature release version, `cglin_scheme_version` is explicitly
`unknown`. Treat comparisons between independently retrieved unknown-version
exports cautiously; use a common frozen export for the reported catalogue.

## Freeze and offline import

`download_cglin_export` writes validated raw batch CSVs and
`cglin_export.json` with requested/returned IDs, explicit missing IDs,
retrieval time, analysis job, content hashes, bytes and attempt counts.
Only fully validated batches receive final filenames. Unexpected IDs,
duplicate output IDs, invalid columns or bad data fail the operation.
Authentication stays in protected configuration or memory and is never
written to exports/manifests; redirect requests receive no API key.

`load_cglin_export(path, scheme=..., scheme_version=..., retrieved_at=...)`
supports UTF-8 CSV, TSV or JSON records (or an `assignments` envelope). It
preserves all original fields in `cglin_export_row`, including raw code,
cgST, available provenance and exact export SHA256. Header-bearing CSV/TSV
must supply an ID column and a code column; sample names cannot replace
source IDs. `annotate_catalogue(rows, assignments)` performs a left join
including every catalogue row and every deduplicated row's source-ID alias.
Identical duplicate imported IDs are counted; differing assignments become
an explicit conflict with no group key. Missing assignments are accounted
for at every depth, rather than collected into a fictitious lineage group.

Flat fields include `cglin_raw`, `cgst`, `cglin_status`,
`cglin_code_status`, `cglin_resolved_depth`, `cglin_provisional`,
`cglin_scheme`, `cglin_scheme_version`, `cglin_retrieved_at`,
`cglin_export_sha256`, `cglin_export_record_count` and
`cglin_group_5/6/7` with corresponding `cglin_status_5/6/7`.

## Focal assignments

`resolve_focal_assignments(focal_rows, annotated_catalogue, crosswalk)`
accepts explicit `{sample_id, source_genome_id}` crosswalk records or
matches exact BioSample, run and assembly accessions. Study IDs, plain
sample names and same-ST membership are not biological identity evidence.
More than one public source match is ambiguous; contradictions between an
explicit crosswalk and accession evidence are conflicts. The returned audit
records evidence, candidates and the outcome for every focal sample.

A supplied focal code requires explicit scheme and provenance (version,
retrieval time or export hash). Valid supplied codes remain unchanged,
including partial ones, with conflicting public evidence reported separately
in `cglin_join_status`. A raw focal code without scheme/provenance receives
`missing_provenance` and no comparable group. Focal samples with no available
assignment receive explicit missing/unavailable/ambiguous status. They are
never uploaded automatically and never borrow a nearest public genome's
complete code. For a novel assembly, ChronoClade can run Pathogenwatch's own
cgMLST caller and `plincer` locally with prepared, pinned databases, or import
validated Pathogenwatch results. BIGSdb-Pasteur submission is required for
definitive identifiers according to its nomenclature guidance. See
[new assembly typing](query-typing.md) for installation and context refinement.

Country figures can still use annotated public records when focal typing is
unavailable. cgLIN annotations do not replace phylogenetic divergence,
recombination or temporal-signal checks.

## Verification

Offline tests cover actual CSV columns and comma codes, complete,
provisional, incomplete, malformed and missing assignments; conflicting
and identical duplicate IDs; different scheme versions; alias-ID joins;
strong/ambiguous/conflicting focal matches; raw export hashes; partial
download batches; corrupt/duplicate/unrequested IDs; retries and credential
redaction. These tests require no live credential.

The 9 October 2026 frozen public same-ST catalogue contained 7,807 raw
records. All requested UUIDs were returned across 79 sequential export
batches in 101.13 seconds, totalling 1,823,172 decoded CSV bytes. Export
manifest SHA256 was
`f6bc0cddd85ba29935e508b9a5bc1dd32dcab62a6804d0b0706d119bdfed9241`.
Its metadata catalogue snapshot hash was
`bcffd21c69089623103170b0e770f6852db3178be0a07a2d029176f94eb4b8d7`.

Raw-record coverage before QC/deduplication was:

| Prefix depth | Groups | Definitive prefix | Provisional prefix | Partial below depth | Missing code |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 | 134 | 7,627 | 152 | 16 | 12 |
| 6 | 285 | 7,627 | 136 | 32 | 12 |
| 7 | 703 | 7,627 | 136 | 32 | 12 |

Every row reconciles to 7,807. No code was malformed. There were 168
nonmissing provisional codes; 12 further records had a provisional cgST
hash but no LIN code. Reimporting the saved exports and repeating annotation
offline produced identical keys. These are **raw genome-record** statistics,
not the QC-filtered, deduplicated geographical counting denominator and not
a claim about validated focal assemblies. The full pilot must compute that
eligible sample-unit denominator and separate focal coverage explicitly.
