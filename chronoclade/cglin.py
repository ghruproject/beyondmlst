"""Frozen Klebsiella cgLIN assignments and conservative focal accession joins.

Depths count components from one; a group includes the *entire* numeric prefix.
No phylogenetic distance or nearest-neighbour code inference occurs here.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener, urlopen

DEFAULT_SCHEME = "scgMLST629_S"
DEFAULT_DEPTHS = (5, 6, 7)
COMPONENT_COUNT = 10


class CGLINError(ValueError):
    """An export or crosswalk cannot safely be interpreted."""


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _field(row: Mapping[str, Any], *names: str) -> Any:
    lookup = {re.sub(r"[^a-z0-9]", "", str(key).lower()): value for key, value in row.items()}
    for name in names:
        value = lookup.get(re.sub(r"[^a-z0-9]", "", name.lower()))
        if value is not None and value != "":
            return value
    return ""


def parse_code(raw: Any) -> tuple[tuple[int, ...], str]:
    """Return the actually resolved prefix and complete/partial/missing/malformed.

    A suffix of unknown components is permitted. Numeric values after an unknown
    component are inconsistent and rejected, rather than spliced into a prefix.
    """
    value = _text(raw)
    if value.lower() in {"", "na", "n/a", "none", "null", "unknown", "-"}:
        return (), "missing"
    tokens = list(raw) if isinstance(raw, (list, tuple)) else re.split(r"[_,]", value)
    if len(tokens) > COMPONENT_COUNT:
        return (), "malformed"
    components: list[int] = []
    unresolved = False
    for token in tokens:
        token = _text(token)
        if re.fullmatch(r"[0-9]+", token):
            if unresolved:
                return (), "malformed"
            components.append(int(token))
        elif token in {"", "-", "?", "*"} or re.fullmatch(r"\*[A-Za-z0-9]+", token):
            unresolved = True
        else:
            return (), "malformed"
    return tuple(components), "complete" if len(components) == COMPONENT_COUNT else "partial"


def group_key(scheme: str, version: str, components: Iterable[int], depth: int) -> str:
    """A stable scheme/version scoped full prefix, never a final component ID."""
    prefix = tuple(components)
    if depth < 1 or depth > COMPONENT_COUNT:
        raise CGLINError("cgLIN depth must be between 1 and 10")
    if len(prefix) < depth:
        return ""
    return json.dumps([scheme, version, list(prefix[:depth])], separators=(",", ":"))


def normalise_assignment(
    row: Mapping[str, Any],
    *,
    scheme: str = DEFAULT_SCHEME,
    scheme_version: str = "unknown",
    retrieved_at: str = "",
    export_sha256: str = "",
    depths: Iterable[int] = DEFAULT_DEPTHS,
) -> dict[str, Any]:
    raw = _field(row, "cglin_raw", "LIN code", "LINcode", "cgLIN", "code")
    cgst = _text(_field(row, "cgst", "cgST", "scgST"))
    components, status = parse_code(raw)
    provisional_value = _field(row, "cglin_provisional", "provisional")
    provisional = (
        provisional_value is True
        or _text(provisional_value).lower() in {"true", "1", "yes"}
        or cgst.startswith("*")
        or "*" in _text(raw)
    )
    scheme = _text(_field(row, "cglin_scheme", "scheme")) or scheme
    version = (
        _text(_field(row, "cglin_scheme_version", "scheme_version", "version")) or scheme_version
    )
    if scheme != DEFAULT_SCHEME:
        components, status = (), "unsupported"
    result = {
        "source_genome_id": _text(_field(row, "source_genome_id", "genome_id", "genome id", "id")),
        "cglin_raw": raw,
        "cgst": cgst,
        "cglin_status": "provisional"
        if provisional and status in {"complete", "partial"}
        else status,
        "cglin_code_status": status,
        "cglin_resolved_depth": len(components),
        "cglin_provisional": provisional,
        "cglin_scheme": scheme,
        "cglin_scheme_version": version,
        "cglin_retrieved_at": _text(row.get("cglin_retrieved_at")) or retrieved_at,
        "cglin_export_sha256": _text(row.get("cglin_export_sha256")) or export_sha256,
        "cglin_export_row": dict(row),
    }
    for depth in depths:
        key = group_key(scheme, version, components, depth)
        result[f"cglin_group_{depth}"] = key
        result[f"cglin_status_{depth}"] = (
            ("provisional" if provisional else "resolved") if key else status
        )
    return result


def load_cglin_export(
    path: str | Path,
    *,
    scheme: str = DEFAULT_SCHEME,
    scheme_version: str = "unknown",
    retrieved_at: str = "",
    depths: Iterable[int] = DEFAULT_DEPTHS,
) -> list[dict[str, Any]]:
    """Import UTF-8 JSON records or header-bearing CSV/TSV; retain exact-byte hash.

    Every record requires an explicit Pathogenwatch ID. Duplicate IDs are kept
    for conflict detection in annotate_catalogue, rather than last-write-wins.
    """
    data = Path(path).read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    text = data.decode("utf-8-sig")
    if text.lstrip().startswith(("[", "{")):
        payload = json.loads(text)
        rows = (
            payload.get("assignments", payload.get("rows"))
            if isinstance(payload, dict)
            else payload
        )
        if not isinstance(rows, list):
            raise CGLINError("cgLIN JSON export must contain a list of records")
    else:
        first_line = text.splitlines()[0] if text.splitlines() else ""
        reader = csv.DictReader(io.StringIO(text), delimiter="\t" if "\t" in first_line else ",")
        headers = {re.sub(r"[^a-z0-9]", "", str(key).lower()) for key in reader.fieldnames or []}
        if not headers.intersection(
            {"sourcegenomeid", "genomeid", "id"}
        ) or not headers.intersection({"cglinraw", "lincode", "cglin", "code"}):
            raise CGLINError("cgLIN CSV/TSV requires source genome ID and LIN code columns")
        rows = list(reader)
        if not first_line:
            raise CGLINError("cgLIN export is empty")
    assignments = []
    for index, row in enumerate(rows, 1):
        if not isinstance(row, dict):
            raise CGLINError(f"cgLIN record {index} is not an object")
        assignment = normalise_assignment(
            row,
            scheme=scheme,
            scheme_version=scheme_version,
            retrieved_at=retrieved_at,
            export_sha256=digest,
            depths=depths,
        )
        if not assignment["source_genome_id"]:
            raise CGLINError(f"cgLIN record {index} has no source genome ID; name joins are unsafe")
        assignments.append(assignment)
    return assignments


def _signature(row: Mapping[str, Any]) -> tuple[str, ...]:
    components, status = parse_code(row.get("cglin_raw"))
    return (
        str(components),
        status,
        *tuple(
            _text(row.get(key))
            for key in ("cgst", "cglin_provisional", "cglin_scheme", "cglin_scheme_version")
        ),
    )


def annotate_catalogue(
    rows: Iterable[Mapping[str, Any]],
    assignments: Iterable[Mapping[str, Any]],
    depths: Iterable[int] = DEFAULT_DEPTHS,
) -> list[dict[str, Any]]:
    """Left join all catalogue records; conflicting exports leave no group key."""
    depths = tuple(depths)
    by_id: dict[str, list[dict[str, Any]]] = {}
    for row in assignments:
        assignment = (
            dict(row) if "cglin_status" in row else normalise_assignment(row, depths=depths)
        )
        ident = _text(assignment.get("source_genome_id"))
        if not ident:
            raise CGLINError("Assignment has no source genome ID")
        by_id.setdefault(ident, []).append(assignment)
    output = []
    for row in rows:
        result = dict(row)
        ident = _text(row.get("source_genome_id"))
        ids = row.get("source_genome_ids") or [ident]
        matches = [a for item in ids for a in by_id.get(_text(item), [])]
        annotation = (
            dict(matches[0])
            if matches
            else normalise_assignment({"source_genome_id": ident}, depths=depths)
        )
        annotation.pop("source_genome_id", None)
        annotation["cglin_export_record_count"] = len(matches)
        if len({_signature(a) for a in matches}) > 1:
            annotation.update(
                cglin_status="conflict",
                cglin_code_status="conflict",
                cglin_resolved_depth=0,
                cglin_conflicting_assignments=matches,
            )
            for depth in depths:
                annotation[f"cglin_group_{depth}"] = ""
                annotation[f"cglin_status_{depth}"] = "conflict"
        result.update(annotation)
        output.append(result)
    return output


_STRONG_ACCESSION = re.compile(r"(?:SAM[NED][A-Z]?\d+|[SED]RR\d+|GC[AF]_\d+(?:\.\d+)?)", re.I)


def strong_accessions(row: Mapping[str, Any]) -> set[str]:
    """Use BioSample/run/assembly evidence, excluding study IDs and plain names."""
    found: set[str] = set()
    for key in (
        "biosample",
        "biosample_accession",
        "run_accessions",
        "assembly_accessions",
        "accession",
        "sample_id",
        "aliases",
    ):
        values = row.get(key) or []
        if not isinstance(values, (list, tuple, set)):
            values = re.split(r"[;,\s]+", str(values))
        for value in values:
            value = _text(value).upper()
            if _STRONG_ACCESSION.fullmatch(value):
                found.add(value)
    return found


def resolve_focal_assignments(
    focal_rows: Iterable[Mapping[str, Any]],
    annotated_catalogue: Iterable[Mapping[str, Any]],
    crosswalk: Iterable[Mapping[str, Any]] | None = None,
    *,
    depths: Iterable[int] = DEFAULT_DEPTHS,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Resolve focal rows through explicit IDs or strong accessions, auditing every join.

    Explicit crosswalk records require sample_id and source_genome_id. More than
    one matching public ID is ambiguous even if its codes happen to agree.
    Supplied focal codes are retained and conflicts explicitly flagged.
    """
    catalogue = list(annotated_catalogue)
    depths = tuple(depths)
    ids: dict[str, Mapping[str, Any]] = {}
    accession_index: dict[str, set[str]] = {}
    for row in catalogue:
        ident = _text(row.get("source_genome_id"))
        if ident in ids:
            raise CGLINError("Focal join requires catalogue IDs to be unique")
        ids[ident] = row
        for accession in strong_accessions(row):
            accession_index.setdefault(accession, set()).add(ident)
    explicit: dict[str, set[str]] = {}
    for link in crosswalk or []:
        sample = _text(link.get("sample_id"))
        ident = _text(link.get("source_genome_id"))
        if not sample or not ident:
            raise CGLINError("Crosswalk requires sample_id and source_genome_id")
        explicit.setdefault(sample, set()).add(ident)
    output, audit = [], []
    seen: set[str] = set()
    for focal in focal_rows:
        sample = _text(focal.get("sample_id"))
        if not sample or sample in seen:
            raise CGLINError("Focal sample IDs must be nonempty and unique")
        seen.add(sample)
        evidence = strong_accessions(focal)
        accession_matches = set().union(*(accession_index.get(a, set()) for a in evidence))
        matches = explicit.get(sample, set()) or accession_matches
        direct_id = _text(focal.get("source_genome_id"))
        if direct_id:
            matches = matches | {direct_id}
        status = "missing" if not matches else "matched" if len(matches) == 1 else "ambiguous"
        if explicit.get(sample) and accession_matches and explicit[sample] != accession_matches:
            status = "conflict"
        if matches - ids.keys() and status == "matched":
            status = "unavailable"
        own_raw = _field(focal, "cglin_raw", "LINcode", "cgLIN")
        own = normalise_assignment(focal, depths=depths)
        provenance = bool(_field(focal, "cglin_scheme", "scheme")) and bool(
            _field(
                focal,
                "cglin_export_sha256",
                "cglin_retrieved_at",
                "cglin_scheme_version",
                "scheme_version",
            )
        )
        own_valid = (
            own["cglin_code_status"] in {"complete", "partial"} and bool(own_raw) and provenance
        )
        if own_raw and not provenance:
            own["cglin_status"] = "missing_provenance"
            for depth in depths:
                own[f"cglin_group_{depth}"] = ""
                own[f"cglin_status_{depth}"] = "missing_provenance"
        matched = ids[next(iter(matches))] if status == "matched" else None
        if matched and own_valid:
            own_prefix, _ = parse_code(own["cglin_raw"])
            matched_prefix, _ = parse_code(matched.get("cglin_raw"))
            overlap = min(len(own_prefix), len(matched_prefix))
            if (
                own["cglin_scheme"] != matched.get("cglin_scheme")
                or own["cglin_scheme_version"] != matched.get("cglin_scheme_version")
                or own_prefix[:overlap] != matched_prefix[:overlap]
                or (own["cgst"] and matched.get("cgst") and own["cgst"] != matched["cgst"])
            ):
                status = "conflict"
        result = dict(focal)
        if own_raw or not matched:
            result.update({k: v for k, v in own.items() if k != "source_genome_id"})
        else:
            result.update(
                {k: v for k, v in matched.items() if k.startswith("cglin_") or k == "cgst"}
            )
        result["cglin_join_status"] = status
        result["cglin_matched_source_genome_id"] = (
            next(iter(matches)) if status == "matched" else ""
        )
        if status in {"ambiguous", "conflict"} and not own_valid:
            result["cglin_status"] = status
            for depth in depths:
                result[f"cglin_group_{depth}"] = ""
                result[f"cglin_status_{depth}"] = status
        output.append(result)
        audit.append(
            {
                "sample_id": sample,
                "status": status,
                "accessions": sorted(evidence),
                "candidate_source_genome_ids": sorted(matches),
                "retained_focal_assignment": own_valid,
            }
        )
    return output, audit


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, *args: Any, **kwargs: Any) -> None:
        return None


def download_cglin_export(
    rows: Iterable[Mapping[str, Any]],
    destination: str | Path,
    *,
    api_key: str,
    job: str = "lincodes-3390273-2",
    scheme_version: str = "unknown",
    base_url: str = "https://pathogen.watch",
    batch_size: int = 100,
    retries: int = 2,
    timeout: float = 60,
    transport: Any = None,
    sleep: Any = time.sleep,
) -> dict[str, Any]:
    """Download bounded analysis exports using the verified numeric-ID contract.

    The API `job` is an existing analysis name, not an asynchronous export ID.
    Store raw decoded CSV batches and a hash-bearing manifest. Join returned
    UUIDs, detect duplicates/unrequested IDs, and account for omitted IDs.
    Authentication is sent only to the API, never its download redirect.
    An injectable transport receives numeric IDs and job and returns CSV bytes.
    """
    rows = list(rows)
    if not api_key:
        raise CGLINError("Pathogenwatch cgLIN export requires an API key")
    if batch_size < 1 or batch_size > 100 or retries < 0:
        raise CGLINError("cgLIN batch size must be 1–100 and retries must be nonnegative")
    if not re.fullmatch(r"lincodes-[A-Za-z0-9-]+", job):
        raise CGLINError("cgLIN job must be an explicit supported lincodes analysis name")
    numeric_to_uuid: dict[str, str] = {}
    for row in rows:
        ident = _text(row.get("source_genome_id"))
        numeric = _text(row.get("numeric_source_id"))
        if not ident or not numeric.isdigit():
            raise CGLINError("cgLIN download requires source_genome_id and numeric_source_id")
        if numeric in numeric_to_uuid or ident in numeric_to_uuid.values():
            raise CGLINError("cgLIN download input IDs must be unique")
        numeric_to_uuid[numeric] = ident
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    assignments, batches = [], []
    retrieved_at = datetime.now(timezone.utc).isoformat()
    items = sorted(numeric_to_uuid, key=int)
    for start in range(0, len(items), batch_size):
        numeric_ids = items[start : start + batch_size]
        requested = {numeric_to_uuid[item] for item in numeric_ids}
        for attempt in range(retries + 1):
            try:
                if transport is not None:
                    data = transport(numeric_ids, job)
                else:
                    query = urlencode({"job": job})
                    request = Request(
                        base_url.rstrip("/") + "/api/downloads/klebsiella-lincodes?" + query,
                        data=json.dumps({"ids": ",".join(numeric_ids)}).encode(),
                        headers={"Content-Type": "application/json", "X-API-Key": api_key},
                        method="POST",
                    )
                    try:
                        response = build_opener(_NoRedirect).open(request, timeout=timeout)
                    except HTTPError as redirect:
                        if redirect.code != 302:
                            raise
                        location = urljoin(base_url, redirect.headers.get("Location", ""))
                        host = (urlsplit(location).hostname or "").lower()
                        # Live verified endpoint uses S3. Reject unexpected redirect
                        # targets; no credential-bearing automatic redirects.
                        if urlsplit(location).scheme != "https" or not (
                            host == urlsplit(base_url).hostname or host.endswith(".amazonaws.com")
                        ):
                            raise CGLINError("cgLIN export returned an unsupported redirect host")
                        response = urlopen(
                            Request(location, headers={"Accept": "text/csv"}), timeout=timeout
                        )
                    with response:
                        data = response.read()
                if data.startswith(b"\x1f\x8b"):
                    data = gzip.decompress(data)
                break
            except HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504} or attempt == retries:
                    raise CGLINError(f"Pathogenwatch cgLIN export HTTP {error.code}") from None
            except (URLError, TimeoutError, OSError):
                if attempt == retries:
                    raise CGLINError("Pathogenwatch cgLIN export transport failed") from None
            sleep(min(2**attempt, 30))
        batch_path = destination / f"batch_{start // batch_size:04d}.csv"
        # Only a fully validated batch becomes a frozen completed file.
        temporary = batch_path.with_suffix(".csv.tmp")
        temporary.write_bytes(data)
        try:
            parsed = load_cglin_export(
                temporary, scheme_version=scheme_version, retrieved_at=retrieved_at
            )
            exported = [record["source_genome_id"] for record in parsed]
            if len(exported) != len(set(exported)):
                raise CGLINError("cgLIN export contains duplicate source IDs")
            if set(exported) - requested:
                raise CGLINError("cgLIN export contains unrequested source IDs")
            temporary.replace(batch_path)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
        assignments.extend(parsed)
        batches.append(
            {
                "file": batch_path.name,
                "sha256": hashlib.sha256(data).hexdigest(),
                "numeric_request_ids": numeric_ids,
                "requested_source_genome_ids": sorted(requested),
                "returned_source_genome_ids": exported,
                "missing_source_genome_ids": sorted(requested - set(exported)),
                "attempts": attempt + 1,
                "bytes": len(data),
            }
        )
    returned = {record["source_genome_id"] for record in assignments}
    manifest = {
        "schema_version": 1,
        "source": "pathogenwatch",
        "job": job,
        "scheme": DEFAULT_SCHEME,
        "scheme_version": scheme_version,
        "retrieved_at": retrieved_at,
        "requested_count": len(rows),
        "returned_count": len(returned),
        "batches": batches,
        "missing_source_genome_ids": sorted(set(numeric_to_uuid.values()) - returned),
        "assignments": assignments,
        "complete": len(returned) == len(rows),
    }
    target = destination / "cglin_export.json"
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, sort_keys=True, indent=2), encoding="utf-8")
    temporary.replace(target)
    return manifest
