"""Bounded, authenticated assembly acquisition with verified content caches."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import shutil
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener


class DownloadError(RuntimeError):
    """A download failed validation or exhausted its retries."""


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def request_download(
    route: str,
    *,
    api_key: str,
    body: dict | None = None,
    base_url: str = "https://pathogen.watch",
    attempts: int = 3,
    timeout: float = 120,
    stats: dict | None = None,
) -> bytes:
    """Follow the documented file redirect without forwarding authentication."""
    if stats is None:
        stats = {}
    for attempt in range(attempts):
        try:
            stats["requests"] = stats.get("requests", 0) + 1
            request = Request(
                base_url + route,
                data=json.dumps(body).encode() if body else None,
                headers={"X-API-Key": api_key, "Content-Type": "application/json"},
            )
            try:
                with build_opener(_NoRedirect()).open(request, timeout=timeout) as response:
                    return response.read()
            except HTTPError as error:
                if error.code not in {301, 302, 303, 307, 308}:
                    raise
                target = urljoin(base_url, error.headers.get("Location", ""))
                parsed = urlparse(target)
                if parsed.scheme != "https" or not parsed.hostname or parsed.username:
                    raise DownloadError("Download redirect is not a valid HTTPS file URL")
                # Signed URLs must never enter logs or manifests. No auth on file request.
                stats["requests"] = stats.get("requests", 0) + 1
                with build_opener().open(Request(target), timeout=timeout) as response:
                    return response.read()
        except HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504}:
                raise DownloadError(f"Pathogenwatch download HTTP {error.code}") from None
        except (URLError, TimeoutError, OSError):
            pass
        if attempt + 1 < attempts:
            time.sleep(min(2**attempt, 8))
    raise DownloadError("Pathogenwatch download exhausted retries")


def validate_fasta(content: bytes) -> bytes:
    """Reject empty, malformed and non-nucleotide assembly files."""
    if content.startswith(b"\x1f\x8b"):
        try:
            content = gzip.decompress(content)
        except (OSError, EOFError) as error:
            raise DownloadError("Corrupt compressed FASTA") from error
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise DownloadError("FASTA is not ASCII") from error
    if not lines or not lines[0].startswith(">"):
        raise DownloadError("Downloaded content is not FASTA")
    length = 0
    current = 0
    seen_header = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header and not current:
                raise DownloadError("FASTA contains an empty contig")
            current = 0
            seen_header = True
        elif re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
            current += len(line)
            length += len(line)
        else:
            raise DownloadError("FASTA contains invalid sequence characters")
    if not length or not current:
        raise DownloadError("FASTA contains no sequence")
    return content


def unpack_bulk(content: bytes, id_to_source: dict[str, str]) -> dict[str, bytes]:
    """Map archive members only by exact IDs, never archive order or fuzzy names."""
    if content.startswith(b"\x1f\x8b"):
        content = gzip.decompress(content)
    try:
        archive = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile as error:
        raise DownloadError("Bulk download is not a ZIP archive") from error
    result: dict[str, bytes] = {}
    for member in archive.infolist():
        if member.is_dir():
            continue
        name = Path(member.filename).name
        stem = re.sub(r"\.(?:fasta|fna|fa)(?:\.gz)?$", "", name, flags=re.I)
        source = id_to_source.get(stem)
        if source is None:
            raise DownloadError("Bulk archive member cannot be joined to a requested ID")
        if source in result:
            raise DownloadError("Bulk archive contains a duplicate genome ID")
        if member.file_size > 100_000_000:
            raise DownloadError("Bulk archive member exceeds the assembly size limit")
        result[source] = validate_fasta(archive.read(member))
    return result


def _atomic(path: Path, content: bytes) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".partial")
    temporary.write_bytes(content)
    temporary.replace(path)


def download_assemblies(
    rows: list[dict],
    *,
    output: Path,
    cache_dir: Path,
    api_key: str,
    workers: int = 2,
    base_url: str = "https://pathogen.watch",
    fetch=None,
) -> tuple[dict[str, Path], list[dict]]:
    """Use verified single-ID downloads; resume only checksum-valid cached files."""
    output.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    def one(row):
        source = str(row["source_genome_id"])
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", source):
            return (
                source,
                None,
                {"source_genome_id": source, "status": "failed", "reason": "unsafe_source_id"},
            )
        metadata = cache_dir / f"{source}.json"
        started = time.monotonic()
        cached = False
        requests = 0
        try:
            data = None
            if metadata.is_file():
                try:
                    saved = json.loads(metadata.read_text())
                    blob = cache_dir / (saved["sha256"] + ".fasta")
                    existing = blob.read_bytes()
                    if (
                        saved["source_genome_id"] == source
                        and saved.get("source_checksum") == row.get("source_checksum")
                        and hashlib.sha256(existing).hexdigest() == saved["sha256"]
                    ):
                        data = validate_fasta(existing)
                        cached = True
                except (ValueError, KeyError, OSError, DownloadError):
                    pass
            if data is None:
                stats = {}
                if fetch:
                    data = validate_fasta(fetch(source))
                    requests = 1
                else:
                    data = validate_fasta(
                        request_download(
                            f"/api/genomes/download/{source}",
                            api_key=api_key,
                            base_url=base_url,
                            stats=stats,
                        )
                    )
                    requests = stats.get("requests", 0)
            digest = hashlib.sha256(data).hexdigest()
            blob = cache_dir / f"{digest}.fasta"
            if not cached:
                _atomic(blob, data)
                _atomic(
                    metadata,
                    json.dumps(
                        {
                            "source_genome_id": source,
                            "sha256": digest,
                            "source_checksum": row.get("source_checksum"),
                        }
                    ).encode(),
                )
            destination = output / f"{source}.fasta"
            shutil.copyfile(blob, destination.with_suffix(".fasta.partial"))
            destination.with_suffix(".fasta.partial").replace(destination)
            return (
                source,
                destination.resolve(),
                {
                    "source_genome_id": source,
                    "status": "cached" if cached else "downloaded",
                    "sha256": digest,
                    "bytes": len(data),
                    "requests": requests,
                    "seconds": round(time.monotonic() - started, 3),
                },
            )
        except (DownloadError, OSError):
            return (
                source,
                None,
                {
                    "source_genome_id": source,
                    "status": "failed",
                    "reason": "download_or_validation_failed",
                    "requests": requests,
                    "seconds": round(time.monotonic() - started, 3),
                },
            )

    with ThreadPoolExecutor(max_workers=max(1, min(workers, 4))) as executor:
        results = list(executor.map(one, rows))
    return (
        {source: path for source, path, _ in results if path is not None},
        [ledger for _, _, ledger in results],
    )
