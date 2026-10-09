import json
from pathlib import Path

import pytest

from chronoclade import typing_setup


def test_species_aliases_and_validation():
    assert typing_setup._normalise_species("K. pneumoniae") == ["klebsiella"]
    assert typing_setup._normalise_species(["ecoli", "Escherichia coli"]) == ["ecoli"]
    assert typing_setup._normalise_species("both") == ["klebsiella", "ecoli"]
    with pytest.raises(typing_setup.TypingSetupError, match="Unsupported"):
        typing_setup._normalise_species("salmonella")


def test_setup_writes_truthful_not_ready_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda name: "/usr/bin/uv")
    monkeypatch.setattr(
        typing_setup,
        "_install_mlst",
        lambda root: {
            "status": "installed",
            "version": typing_setup.UPSTREAM["mlst"]["version"],
            "commit": typing_setup.UPSTREAM["mlst"]["commit"],
            "command": ["node", "index.js"],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_install_python_tool",
        lambda uv, root, name, python: {
            "status": "installed",
            "version": typing_setup.UPSTREAM[name]["version"],
            "commit": typing_setup.UPSTREAM[name]["commit"],
            "command": [name],
        },
    )

    manifest_path, report = typing_setup.setup_typing(tmp_path, species="ecoli")

    assert manifest_path == tmp_path / "typing_setup.json"
    assert json.loads(manifest_path.read_text()) == report
    assert report["species"] == ["ecoli"]
    assert report["tools"]["hclink"]["status"] == "installed"
    assert report["run_ready"] is False
    assert report["query_config_path"] == str(tmp_path / "query_config.json")
    query_config = json.loads(Path(report["query_config_path"]).read_text())
    assert query_config["ready"] is False
    assert query_config["organisms"]["Escherichia coli"]["cgmlst"]["scheme_version"] is None
    assert query_config["organisms"]["Escherichia coli"]["assignment"]["loci_file"].endswith(
        "/databases/hclink/loci.json"
    )
    db = report["databases"]["ecoli"]
    assert db["cgmlst_ready"] is False
    assert db["lineage_ready"] is False
    assert "2024-12-31" in db["public_cgmlst_data"]
    assert "api key" in db["lineage_database_authentication"].lower()


def test_prepare_db_is_explicit_and_records_builds(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda name: "/usr/bin/uv")
    monkeypatch.setattr(
        typing_setup,
        "_install_mlst",
        lambda root: {
            "status": "installed",
            "version": typing_setup.UPSTREAM["mlst"]["version"],
            "commit": typing_setup.UPSTREAM["mlst"]["commit"],
            "command": ["node", "index.js"],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_install_python_tool",
        lambda uv, root, name, python: {
            "status": "installed",
            "version": typing_setup.UPSTREAM[name]["version"],
            "commit": typing_setup.UPSTREAM[name]["commit"],
            "command": [name],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_prepare_public_cgmlst",
        lambda uv, root, species, tools: {
            "status": "built",
            "scheme": typing_setup.SPECIES[species]["mlst_scheme"],
        },
    )

    _, report = typing_setup.setup_typing(tmp_path, species="klebsiella", prepare_db=True)

    assert report["cgmlst_builds"]["klebsiella"] == {
        "status": "built",
        "scheme": "klebsiella_1",
    }
    assert "only when their credential files were supplied" in report["database_preparation"]


def test_missing_uv_is_reported_without_writing(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda _name: None)
    with pytest.raises(typing_setup.TypingSetupError, match="uv executable is unavailable"):
        typing_setup.setup_typing(tmp_path)
    assert not (tmp_path / "typing_setup.json").exists()


def test_hclink_builder_reads_key_from_file_without_putting_it_in_command(tmp_path, monkeypatch):
    key_file = tmp_path / "key.txt"
    key_file.write_text("private-fixture-key")
    tool = {"source": str(tmp_path / "source"), "environment": str(tmp_path / "venv")}

    def run(args, **kwargs):
        assert all("private-fixture-key" not in str(arg) for arg in args)
        assert kwargs["env"]["CHRONOCLADE_ENTEROBASE_API_KEY_FILE"] == str(key_file)
        script = Path(args[1]).read_text()
        assert "private-fixture-key" not in script
        db = Path(kwargs["env"]["CHRONOCLADE_HCLINK_DB"])
        (db / "metadata.json").write_text(json.dumps({"datestamp": "2026-01-01T12:00:00"}))
        (db / "loci.json").write_text(json.dumps({"genes": ["a", "b"]}))
        (db / "index.usearch").write_bytes(b"controlled-fixture-index")
        (db / "alleles.db").write_bytes(b"fixture-alleles")
        (db / "ST.txt.xz").write_bytes(b"fixture-clusters")
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._build_hclink_db(tmp_path, tool, key_file)
    assert report["scheme_version"] == "2026-01-01"
    assert "private-fixture-key" not in json.dumps(report)


def test_plincer_credential_cache_is_private_and_outside_reference_database(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "scheme.toml").write_text("controlled-fixture")
    secret = tmp_path / "secrets.json"
    secret.write_text("{}")

    def run(args, **kwargs):
        cache = Path(args[args.index("--secrets-cache-file") + 1])
        assert cache.parent == tmp_path / "credential-cache"
        assert cache.parent.stat().st_mode & 0o777 == 0o700
        db = tmp_path / "databases" / "plincer"
        (db / "profiles.json.xz").write_bytes(b"fixture")
        (db / "alleles.sqlite").write_bytes(b"fixture")
        (db / "metadata.json").write_text(
            json.dumps({"last_updated": "2026-01-01", "genes": ["a"]})
        )
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._build_plincer_db(
        tmp_path, {"source": str(source), "command": ["plincer"]}, secret
    )
    assert report["status"] == "built"
    assert report["scheme_version"] == "2026-01-01"


def test_public_scheme_indexer_receives_its_own_local_scheme_manifest(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    interpreter = tmp_path / "venvs" / "typing_databases" / "bin" / "python"
    interpreter.parent.mkdir(parents=True)
    interpreter.write_text("fixture")
    (interpreter.parent / "download_schemes").write_text("fixture")
    monkeypatch.setattr(typing_setup, "_clone_at", lambda *args: source)
    monkeypatch.setattr(typing_setup, "_python_env", lambda *args: interpreter)
    other = tmp_path / "databases" / "mlst-schemes" / "ecoli_1-selection.json"
    other.parent.mkdir(parents=True)
    other.write_text("other-scheme-preserved")
    calls = []

    def run(args, **kwargs):
        calls.append(args)
        if "--output-schemes-file" in args:
            folder = Path(args[args.index("--output-dir") + 1]) / "pasteur" / "klebsiella_1"
            folder.mkdir(parents=True)
            (folder / "metadata.json").write_text(json.dumps({"genes": ["a"]}))
            selected = Path(args[args.index("--output-schemes-file") + 1])
            selected.write_text(
                json.dumps(
                    {
                        "schemes": [
                            {
                                "shortname": "klebsiella_1",
                                "db_path": str(folder),
                                "last_updated": "2026-01-01",
                            }
                        ]
                    }
                )
            )
        if "index" in args:
            database = Path(
                next(arg.split("=", 1)[1] for arg in args if arg.startswith("--database="))
            )
            manifest = json.loads((database / "schemes.json").read_text())
            record = manifest["schemes"][0]
            assert record["shortname"] == "klebsiella_1"
            assert record["db_path"] == "."
            assert (database / record["db_path"] / "metadata.json").is_file()
            index = Path(next(arg.split("=", 1)[1] for arg in args if arg.startswith("--index=")))
            index.mkdir()
            (index / "klebsiella_1_metadata.json").write_text("fixture-index")
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._prepare_public_cgmlst(
        "uv",
        tmp_path,
        "klebsiella",
        {"mlst": {"npm": "npm", "source": str(source), "command": ["node", "index.js"]}},
    )
    assert report["status"] == "built"
    assert other.read_text() == "other-scheme-preserved"


def test_hclink_readiness_requires_allele_lookup_and_cluster_info(tmp_path):
    database = tmp_path / "databases" / "hclink"
    database.mkdir(parents=True)
    for name in ("metadata.json", "index.usearch", "loci.json"):
        (database / name).write_text("fixture")
    state = typing_setup._database_state(tmp_path, "ecoli")
    assert not state["lineage_ready"]
    assert str(database / "alleles.db") in state["missing_resources"]
    assert str(database / "ST.txt.xz") in state["missing_resources"]
    for name in ("alleles.db", "ST.txt.xz"):
        (database / name).write_text("fixture")
    assert typing_setup._database_state(tmp_path, "ecoli")["lineage_ready"]
