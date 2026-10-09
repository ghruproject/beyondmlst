import io
import json
import zipfile

import pytest

from chronoclade.pathogenwatch_download import (
    DownloadError,
    download_assemblies,
    unpack_bulk,
    validate_fasta,
)


def test_download_resume_verifies_content_and_recovers_corruption(tmp_path):
    calls = []

    def fetch(source):
        calls.append(source)
        return b">contig\nACGTN\n"

    arguments = dict(
        rows=[{"source_genome_id": "UUID1"}],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=fetch,
    )
    paths, ledger = download_assemblies(**arguments)
    assert paths["UUID1"].read_bytes() == b">contig\nACGTN\n"
    assert ledger[0]["status"] == "downloaded"
    _, ledger = download_assemblies(**arguments)
    assert calls == ["UUID1"] and ledger[0]["status"] == "cached"
    meta = json.loads((tmp_path / "cache/UUID1.json").read_text())
    (tmp_path / "cache" / f"{meta['sha256']}.fasta").write_text("corrupt")
    _, ledger = download_assemblies(**arguments)
    assert calls == ["UUID1", "UUID1"] and ledger[0]["status"] == "downloaded"
    assert "secret" not in json.dumps(ledger)


def test_failures_preserve_successful_downloads(tmp_path):
    paths, ledger = download_assemblies(
        [{"source_genome_id": "a"}, {"source_genome_id": "b"}],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=lambda source: b">ok\nACGT\n" if source == "a" else b"<html>error</html>",
    )
    assert set(paths) == {"a"}
    assert [row["status"] for row in ledger] == ["downloaded", "failed"]
    assert not (tmp_path / "out/b.fasta").exists()


def test_bulk_matches_ids_not_archive_order_and_retains_partial_batch():
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as archive:
        archive.writestr("b.fasta", ">b\nACGT\n")
    assert set(unpack_bulk(stream.getvalue(), {"a": "UUIDa", "b": "UUIDb"})) == {"UUIDb"}
    with pytest.raises(DownloadError, match="joined"):
        unpack_bulk(stream.getvalue(), {"a": "UUIDa"})


@pytest.mark.parametrize("data", [b"", b">a\n", b"<html>", b">a\nACGTX\n", b">a\nACGT\n>b\n"])
def test_invalid_fasta_rejected(data):
    with pytest.raises(DownloadError):
        validate_fasta(data)


def test_identical_contents_download_concurrently(tmp_path):
    paths, ledger = download_assemblies(
        [{"source_genome_id": str(i)} for i in range(4)],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        workers=4,
        fetch=lambda source: b">contig\nACGT\n",
    )
    assert len(paths) == 4
    assert all(row["status"] == "downloaded" for row in ledger)
    assert len(list((tmp_path / "cache").glob("*.fasta"))) == 1


def test_provider_revision_invalidates_cached_assembly(tmp_path):
    calls = []

    def fetch(source):
        calls.append(source)
        return b">a\nACGT\n" if len(calls) == 1 else b">a\nTGCA\n"

    row = {"source_genome_id": "a", "source_checksum": "revision1"}
    kwargs = dict(
        output=tmp_path / "out", cache_dir=tmp_path / "cache", api_key="secret", fetch=fetch
    )
    download_assemblies([row], **kwargs)
    row["source_checksum"] = "revision2"
    paths, ledger = download_assemblies([row], **kwargs)
    assert len(calls) == 2 and ledger[0]["status"] == "downloaded"
    assert b"TGCA" in paths["a"].read_bytes()


def test_wrong_source_content_is_rejected(tmp_path):
    import hashlib

    row = {
        "source_genome_id": "a",
        "source_checksum": hashlib.sha1(b">a\nACGT\n").hexdigest(),
        "source_length": 4,
    }
    paths, ledger = download_assemblies(
        [row],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=lambda source: b">a\nTGCA\n",
    )
    assert not paths and ledger[0]["reason"] == "source_checksum_mismatch"
    assert not list((tmp_path / "cache").glob("*.fasta"))
