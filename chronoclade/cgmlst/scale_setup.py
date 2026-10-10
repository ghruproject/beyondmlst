"""Explicit, pinned source installation for the optional RapidNJ backend."""

from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

RAPIDNJ_REVISION = "ed2d36e219d9db16778b941b5054c0fd021b528a"
RAPIDNJ_URL = f"https://github.com/somme89/rapidNJ/archive/{RAPIDNJ_REVISION}.tar.gz"
RAPIDNJ_SOURCE_SHA256 = "57ea662a3589459528beaa96de011bd1b389630484d9e3215a942bdb19518f21"
DEFAULT_RAPIDNJ = Path.home() / ".local/share/chronoclade/scale/bin/rapidnj"


def rapidnj_executable(requested=None):
    """Find an explicit executable, PATH backend or verified persistent install."""
    explicit = requested or os.environ.get("CHRONOCLADE_RAPIDNJ")
    if explicit:
        binary = shutil.which(str(explicit))
        if binary is None:
            raise RuntimeError(f"Requested RapidNJ executable is unavailable: {explicit}")
        return binary
    binary = shutil.which("rapidnj")
    if binary:
        return binary
    if DEFAULT_RAPIDNJ.is_file() and os.access(DEFAULT_RAPIDNJ, os.X_OK):
        receipt_path = DEFAULT_RAPIDNJ.parent / "rapidnj-install.json"
        if receipt_path.is_file():
            receipt = json.loads(receipt_path.read_text())
            digest = hashlib.sha256(DEFAULT_RAPIDNJ.read_bytes()).hexdigest()
            if (
                receipt.get("binary_sha256") == digest
                and receipt.get("source_revision") == RAPIDNJ_REVISION
            ):
                return str(DEFAULT_RAPIDNJ)
    raise RuntimeError(
        "Large cgMLST blocks require RapidNJ. Run chronoclade setup-scale, install rapidnj on PATH, or set CHRONOCLADE_RAPIDNJ. No tree is substituted or sampled."
    )


def install_rapidnj(target=None, *, jobs=2):
    """Compile untouched pinned upstream sources and save a provenance receipt.

    Upstream uses x86 SSE instructions. Apple Silicon therefore needs an x86_64
    binary and the installed Rosetta runtime. Unsupported ARM/Linux platforms
    fail explicitly instead of patching the upstream distance/NJ algorithms.
    """
    target = Path(target or DEFAULT_RAPIDNJ).expanduser().resolve()
    if not isinstance(jobs, int) or jobs < 1:
        raise ValueError("Build jobs must be a positive integer")
    if target.exists():
        receipt_path = target.parent / "rapidnj-install.json"
        if receipt_path.is_file():
            receipt = json.loads(receipt_path.read_text())
            if (
                receipt.get("source_revision") == RAPIDNJ_REVISION
                and receipt.get("binary_sha256") == hashlib.sha256(target.read_bytes()).hexdigest()
            ):
                return target
        raise ValueError(f"Refusing to overwrite an unverified existing RapidNJ binary: {target}")
    machine, system = platform.machine().lower(), platform.system()
    command = ["make", f"-j{jobs}"]
    if system == "Darwin" and machine in {"arm64", "aarch64"}:
        rosetta = subprocess.run(["arch", "-x86_64", "/usr/bin/true"], capture_output=True)
        if rosetta.returncode:
            raise RuntimeError(
                "RapidNJ upstream requires x86 SSE; this Apple Silicon host needs Rosetta installed before setup-scale can build/run it."
            )
        command += ["CC=clang++ -arch x86_64", "LINK=clang++ -arch x86_64"]
    elif machine not in {"x86_64", "amd64", "i386", "i686"}:
        raise RuntimeError(
            f"Pinned upstream RapidNJ SSE source does not support {system}/{machine}; provide a separately validated executable explicitly."
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="chronoclade-rapidnj-") as work:
        archive = Path(work) / "source.tar.gz"
        with urlopen(RAPIDNJ_URL, timeout=60) as response, archive.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != RAPIDNJ_SOURCE_SHA256:
            raise RuntimeError(
                "RapidNJ source archive checksum mismatch; nothing was compiled or installed"
            )
        with tarfile.open(archive) as source:
            for member in source.getmembers():
                destination = (Path(work) / member.name).resolve()
                if (
                    not destination.is_relative_to(Path(work).resolve())
                    or member.issym()
                    or member.islnk()
                ):
                    raise RuntimeError("Unsafe member in RapidNJ source archive")
            source.extractall(work, filter="data")
        directory = Path(work) / f"rapidNJ-{RAPIDNJ_REVISION}"
        built = subprocess.run(command, cwd=directory, text=True, capture_output=True)
        (target.parent / "rapidnj-build.log").write_text(built.stdout + built.stderr)
        if built.returncode:
            raise RuntimeError(
                f"RapidNJ compilation failed; see {target.parent / 'rapidnj-build.log'}"
            )
        binary = directory / "bin/rapidnj"
        checked = subprocess.run([str(binary), "-h"], text=True, capture_output=True)
        if "Rapid neighbour-joining" not in checked.stdout + checked.stderr:
            raise RuntimeError("Built RapidNJ executable did not pass the startup check")
        # A tiny real distance input exercises the compiled NJ parser/algorithm.
        matrix = Path(work) / "check.phy"
        matrix.write_text("3\nA 0 0.1 0.2\nB 0.1 0 0.3\nC 0.2 0.3 0\n")
        checked = subprocess.run(
            [str(binary), str(matrix), "-i", "pd", "-o", "t"], text=True, capture_output=True
        )
        if checked.returncode or not all(name in checked.stdout for name in ("A", "B", "C")):
            raise RuntimeError("Built RapidNJ failed its three-tip distance-tree check")
        shutil.copy2(binary, target)
    receipt = dict(
        source_url=RAPIDNJ_URL,
        source_revision=RAPIDNJ_REVISION,
        source_archive_sha256=RAPIDNJ_SOURCE_SHA256,
        platform=f"{system}/{machine}",
        build_command=command,
        binary_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        startup_and_three_tip_check="passed",
        binary=str(target),
    )
    (target.parent / "rapidnj-install.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return target
