# Frozen DNA inputs for optional ESM2 preparation

`chronoclade.esm2.allele_sequences.prepare_allele_sequences` translates actual
allele DNA and writes a protein FASTA plus an explicit sample/locus mapping for
the existing ESM2 embedding engine. This preparation step is offline and does
not need PyTorch, model weights, a network connection or Pasteur credentials.
It does not obtain DNA from allele IDs, hashes, proteins or neighbouring samples.

## Catalogue contract

Supply a UTF-8 JSON manifest with these fields:

| Field | Required value |
| --- | --- |
| `schema` | `chronoclade.esm2.allele-sequence-catalogue` |
| `schema_version` | `1` |
| `status` | `complete` |
| `scope.species` | Exact nonempty species text used by every prepared sample |
| `scope.scheme_id` | Exact prepared locus catalogue scheme ID |
| `scope.scheme_version` | Exact prepared scheme version |
| `scope.database_version` | Exact prepared database version, or JSON `null` when both are unknown |
| `scope.database_sha256` | Known lowercase 64-character SHA256 matching the prepared catalogue |
| `artifact.path` | Local CSV path relative to the manifest directory |
| `artifact.sha256` | SHA256 of the exact CSV bytes |

The `scope` object must contain exactly those five fields. The CSV must begin
with exactly `locus,allele,dna`, in that order. Each row contains one locus,
one categorical allele identifier and its DNA sequence; no sample IDs occur in
this input. Each locus must belong to the prepared scheme, and each locus/allele
pair must occur only once. Duplicate records are rejected even when their DNA
is identical. Missing-call markers such as `0` cannot be used as allele IDs.
Leading zeroes and novel allele hashes are preserved as text.

Example CSV rows, **synthetic validation controls**:

```csv
locus,allele,dna
a,001,ATGAAATAA
a,002,ATGAAGTAA
b,1,GTGGCTTAG
```

The first two rows retain distinct allele identities and DNA hashes but translate
to the same protein `MK`. The last row translates to `MA` using bacterial CDS
initiation. These controls are not clinical samples or a production database.

Freeze the genuine database identity and CSV checksum before use. A scheme's
locus-list checksum, a collection of allele IDs or this example must not be
presented as the fingerprint of a different allele database. Unknown prepared
database fingerprints are rejected. A manifest can describe a partial sequence
catalogue: `status=complete` means the frozen input was fully published, not that
it contains every possible allele. Missing called alleles remain explicit.

Exactly one profile matrix is supported per call. Split multiple schemes or
database scopes into separately prepared cohorts. Taxonomic aliases are not
guessed: unknown species, mixed species or differences in species text are
errors. CSV references must remain within the manifest directory; absolute paths,
parent-directory traversal and remote retrieval are unsupported.

## Python API

```python
from pathlib import Path
from chronoclade.datasets import load_dataset
from chronoclade.esm2.allele_sequences import prepare_allele_sequences
from chronoclade.esm2 import run_embeddings

prepared = load_dataset(Path("prepared/dataset.json"))
sequences = prepare_allele_sequences(
    prepared,
    Path("frozen-alleles/catalogue.json"),
    Path("sequence-input"),
    genetic_code=11,
    terminal_stop="require",
    invalid_policy="exclude",
    require_complete=False,
    max_length=1022,
)

# Requires the separately installed optional inference environment and local weights.
embeddings = run_embeddings(
    sequences.fasta_path,
    Path("embeddings"),
    model="8M",
    checkpoint=Path("weights/esm2_t6_8M_UR50D.pt"),
    allow_download=False,
)
```

`SequenceMappingResult` provides `manifest_path`, `fasta_path`, `mapping_path`
and the loaded `manifest` dictionary. Pass `mapping_path` together with the
embedding manifest to the existing sample/locus distance API. The output
directory must be new or empty, so an existing published bundle cannot be
silently overwritten.

Only NCBI bacterial genetic code 11 is supported. DNA is normalised to uppercase;
only `A`, `C`, `G` and `T` are accepted. The length must be in frame and the first
codon must be a valid table-11 initiator. Initiation is translated as methionine,
including alternative initiators, with that action recorded for each allele.
Internal stops, gaps, ambiguous bases, invalid starts, out-of-frame sequences and
proteins exceeding `max_length` are excluded with distinct reasons. No sequence
is truncated, repaired, reverse-complemented or imputed.

The default `terminal_stop="require"` requires a terminal stop and removes that
single codon before embedding. `terminal_stop="allow_absent"` explicitly permits
a sequence without a terminal stop; a present terminal stop is still removed.
This option does not accept internal stops or relax frame/initiation validation.
Only called alleles are translated; unreferenced CSV alleles are counted in the
input evidence but do not become protein records.

## Output and audit

The bundle contains:

- `proteins.fasta`: one validated record per mapped locus/allele identity. Record
  IDs are deterministic hashes of the complete scope, locus and categorical
  allele. Identical proteins can occur in distinct records: the embedding reader
  embeds each unique amino-acid sequence once and preserves every record mapping.
- `sample_loci.csv`: exactly `sample_id,locus,record_id`; only validated mapped
  cells are included. Repeated use of one allele references its existing record.
- `sequence_mapping.json`: versioned status, frozen-input paths and hashes,
  prepared sample/profile/locus hashes, translation parameters, output paths and
  hashes, counts, allele audit, and every sample/locus cell.

Every prepared sample is retained in `samples` and `cells`, including samples
absent from the profile matrix. Cells distinguish `mapped`, `missing_profile`,
`missing_call`, `missing_sequence` and `invalid_cds`. Missing cells have no protein
ID and are never represented as a zero protein/vector. Per-allele audit retains
locus/allele identity, record ID, DNA hash, protein hash and translation actions
when available. Synonymous differences remain visible even when their protein
hashes coincide.

`status="complete"` means the mapping bundle was published successfully;
`coverage="partial"` explicitly records unavailable cells. A bundle can have no
mapped proteins and still provide a complete preparation audit; the embedding
engine then has no usable inference input.

`require_complete=True` requires every sample to have a validated sequence at
every prepared locus. `invalid_policy="error"` fails if any called allele has an
invalid CDS. Either strict failure publishes an audit with `status="failed"`,
then raises `SequenceMappingError` with `manifest_path` pointing to that audit.
Compatibility, structural and checksum errors fail before output publication.
Callers should stop inference on the exception and inspect its saved audit.

## Validation evidence

The tests cover synonymous allele identity, alternative bacterial initiators,
malformed CDSs, explicit stop handling, missing calls/sequences/profiles,
species/scheme/version/fingerprint conflicts, checksums and immutable outputs.
Import checks verify this preparation module loads without PyTorch or fair-esm.

An offline check on 10 October 2026 used the retained actual Pasteur DNA JSON
records in `tmp/esm2-benchmark/alleles/`. All 32 alleles selected in
`validation/esm2/fixture-manifest.json` translated exactly to that fixture's
proteins and matched its DNA SHA256s. This validates translation against real
retained sequences; it does not claim a validated clinical-genome mapping or a
matching complete allele-database snapshot. The committed fixture retains
proteins and DNA hashes, so the actual DNA files are necessary to repeat that
specific check.
