"""Supporting job preparation; submission remains explicit."""
from pathlib import Path
from typing import Annotated
import typer


def job_command(stage: str, arguments: Annotated[list[str], typer.Argument()],
                output: Annotated[Path, typer.Option('--out', '-o')] = Path('stage_job'),
                slurm: bool = False, cpus: int = 4, memory_gb: int = 16,
                walltime: str = '04:00:00', account: str | None = None,
                partition: str | None = None, gpus: int = 0,
                executable: str = 'chronoclade'):
    """Write an exact local/SLURM command; use -- before stage flags."""
    from .jobs import JobSpec, SlurmResources, write_job
    try:
        spec = JobSpec(stage, tuple(arguments), Path.cwd(), executable)
        resources = SlurmResources(cpus, memory_gb, walltime, account, partition, gpus) if slurm else None
        path = write_job(spec, output, resources=resources)
    except (ValueError, OSError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(2) from exc
    typer.echo(f'Prepared job: {path}')
