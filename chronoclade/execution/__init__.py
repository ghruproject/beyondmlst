"""Stage-neutral local and SLURM command execution."""
from .jobs import JobSpec, SlurmResources, write_job

__all__ = ["JobSpec", "SlurmResources", "write_job"]
