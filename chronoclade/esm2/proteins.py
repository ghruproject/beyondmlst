"""Validated protein inputs and bounded, padding-aware batches."""

from dataclasses import dataclass
import hashlib
from pathlib import Path


class EmbeddingError(ValueError):
    """An embedding input, dependency or saved artifact is invalid."""


@dataclass(frozen=True)
class Protein:
    id: str
    sequence: str


@dataclass(frozen=True)
class ProteinInputs:
    proteins: tuple[Protein, ...]
    mapping: tuple[dict, ...]
    exclusions: tuple[dict, ...]
    input_sha256: str
    records: int


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_proteins(path: Path, *, max_length: int = 1022,
                  invalid_policy: str = "error") -> ProteinInputs:
    """Upper-case protein FASTA, retaining every source record or its exclusion."""
    if invalid_policy not in {"error", "exclude"}:
        raise EmbeddingError("invalid_policy must be 'error' or 'exclude'")
    if not 1 <= max_length <= 1022:
        raise EmbeddingError("max_length must be between 1 and 1022; truncation is unsupported")
    path = Path(path)
    if not path.is_file():
        raise EmbeddingError(f"Protein FASTA does not exist: {path}")
    records = []
    header = None
    pieces = []
    try:
        with path.open(encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        records.append((header, "".join(pieces).upper()))
                    header, pieces = line[1:].strip(), []
                    if not header:
                        raise EmbeddingError(f"Empty FASTA header on line {number}")
                elif header is None:
                    raise EmbeddingError(f"Sequence before FASTA header on line {number}")
                else:
                    pieces.append(line)
    except UnicodeError as exc:
        raise EmbeddingError("Protein FASTA must be UTF-8 text") from exc
    if header is not None:
        records.append((header, "".join(pieces).upper()))
    if not records:
        raise EmbeddingError("Protein FASTA contains no records")
    seen = set()
    unique = {}
    mapping, exclusions = [], []
    residues = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
    for header, sequence in records:
        record_id = header.split()[0]
        if record_id in seen:
            raise EmbeddingError(f"Duplicate FASTA record ID: {record_id}")
        seen.add(record_id)
        reason = None
        if not sequence:
            reason = "empty sequence"
        elif set(sequence) - residues:
            reason = f"unsupported residues: {''.join(sorted(set(sequence) - residues))}"
        elif len(sequence) > max_length:
            reason = f"length {len(sequence)} exceeds max_length {max_length}"
        record = {"record_id": record_id, "header": header, "length": len(sequence)}
        if reason:
            if invalid_policy == "error":
                raise EmbeddingError(f"Invalid protein {record_id}: {reason}; no truncation applied")
            exclusions.append({**record, "reason": reason})
            continue
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        unique.setdefault(digest, Protein(digest, sequence))
        mapping.append({**record, "protein_id": digest})
    if not unique:
        raise EmbeddingError("No valid protein sequences remain after exclusions")
    return ProteinInputs(tuple(unique.values()), tuple(mapping), tuple(exclusions),
                         file_sha256(path), len(records))


def token_batches(proteins, token_budget: int):
    """Bound padded batch token count, including one BOS and one EOS token."""
    if token_budget < 3:
        raise EmbeddingError("token_budget must be at least 3")
    batch = []
    width = 0
    for protein in sorted(proteins, key=lambda p: (len(p.sequence), p.id)):
        tokens = len(protein.sequence) + 2
        if tokens > token_budget:
            raise EmbeddingError(
                f"Protein {protein.id} needs {tokens} tokens, exceeding token_budget {token_budget}"
            )
        next_width = max(width, tokens)
        if batch and next_width * (len(batch) + 1) > token_budget:
            yield tuple(batch)
            batch, width = [], 0
        batch.append(protein)
        width = max(width, tokens)
    if batch:
        yield tuple(batch)
