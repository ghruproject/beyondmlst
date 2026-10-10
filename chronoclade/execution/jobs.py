"""Portable explicit jobs: one scientific command for local and SLURM runs."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re
import shlex
import subprocess

from chronoclade.artifacts import file_sha256, write_json


@dataclass(frozen=True)
class SlurmResources:
    cpus: int = 4
    memory_gb: int = 16
    walltime: str = "04:00:00"
    account: str | None = None
    partition: str | None = None
    gpus: int = 0

    def __post_init__(self):
        if any(type(v) is not int for v in (self.cpus, self.memory_gb, self.gpus)):
            raise ValueError("Resource counts must be integers")
        if self.cpus < 1 or self.memory_gb < 1 or self.gpus < 0:
            raise ValueError("CPU/memory must be positive and GPU count nonnegative")
        if not re.fullmatch(r"(?:\d+-)?\d{1,3}:\d{2}:\d{2}", self.walltime):
            raise ValueError("walltime must use [days-]hours:minutes:seconds")
        if any(int(part) > 59 for part in self.walltime.split(":")[1:]):
            raise ValueError("Invalid walltime minutes/seconds")
        for value in (self.account, self.partition):
            if value is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
                raise ValueError("SLURM account/partition must be a single safe name")


@dataclass(frozen=True)
class JobSpec:
    stage: str
    arguments: tuple[str, ...]
    working_directory: Path
    executable: str = "chronoclade"

    def __post_init__(self):
        if self.stage not in {"prepare", "cgmlst", "esm2", "tree", "time", "time-compare"}:
            raise ValueError("Unknown analysis stage")
        if not self.working_directory.is_dir():
            raise ValueError("Job working directory does not exist")
        if not self.executable or "\x00" in self.executable:
            raise ValueError("Job executable must be nonempty")
        if any(not isinstance(arg, str) or "\x00" in arg for arg in self.arguments):
            raise ValueError("Job arguments must be strings without NUL characters")
        if any(any(word in arg.lower() for word in ("api-key", "api_key", "token", "password"))
               for arg in self.arguments):
            raise ValueError("Credentials must come from runtime configuration, not job arguments")

    @property
    def argv(self):
        return [self.executable, self.stage, *self.arguments]


def write_job(spec: JobSpec, output: Path, *, resources=None):
    """Save script and manifest, without submitting or running the job."""
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Job output must be empty")
    output.mkdir(parents=True, exist_ok=True)
    directives = []
    if resources is not None:
        directives = [f"#SBATCH --cpus-per-task={resources.cpus}",
                      f"#SBATCH --mem={resources.memory_gb}G", f"#SBATCH --time={resources.walltime}"]
        for key in ("account", "partition"):
            if getattr(resources, key):
                directives.append(f"#SBATCH --{key}={getattr(resources, key)}")
        if resources.gpus:
            directives.append(f"#SBATCH --gpus={resources.gpus}")
    script = output / ("run.slurm" if resources else "run.sh")
    script.write_text("#!/usr/bin/env bash\n" + "\n".join(directives) + "\nset -euo pipefail\n"
                      + "cd -- " + shlex.quote(str(spec.working_directory.resolve())) + "\n"
                      + "exec " + shlex.join(spec.argv) + "\n")
    script.chmod(0o700)
    manifest = output / "job.json"
    write_json(manifest, dict(schema="chronoclade.execution.job", schema_version=1,
                             status="prepared", stage=spec.stage, argv=spec.argv,
                             working_directory=str(spec.working_directory.resolve()),
                             executor="slurm" if resources else "local",
                             resources=asdict(resources) if resources else None,
                             script=dict(path=script.name, sha256=file_sha256(script))))
    return manifest


def run_local(spec: JobSpec, *, log_path: Path):
    """Execute exactly the declared argv; preserve the process return code and log."""
    with Path(log_path).open("w") as log:
        process = subprocess.run(spec.argv, cwd=spec.working_directory,
                                 stdout=log, stderr=subprocess.STDOUT, check=False)
    return process.returncode
