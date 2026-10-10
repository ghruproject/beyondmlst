import subprocess
from pathlib import Path
import pytest
from chronoclade.execution.jobs import JobSpec, SlurmResources, write_job


def test_script_quotes_stage_arguments_and_no_submission(tmp_path):
    target = tmp_path / 'job'
    spec = JobSpec('cgmlst', ('data with spaces.json', '--out', '$(touch exploit)'), tmp_path)
    write_job(spec, target, resources=SlurmResources(account='team', partition='cpu'))
    script = (target / 'run.slurm').read_text()
    assert "'$(touch exploit)'" in script
    assert '#SBATCH --account=team' in script and 'sbatch ' not in script
    assert subprocess.run(['bash', '-n', str(target / 'run.slurm')]).returncode == 0
    with pytest.raises(ValueError, match='empty'):
        write_job(spec, target)


@pytest.mark.parametrize('resources', [dict(account='bad\n#SBATCH --foo'), dict(cpus=0),
                                      dict(walltime='01:99:00'), dict(gpus=-1)])
def test_invalid_resources(resources):
    with pytest.raises(ValueError):
        SlurmResources(**resources)


def test_credentials_rejected(tmp_path):
    with pytest.raises(ValueError, match='Credentials'):
        JobSpec('prepare', ('--api-key', 'secret'), Path(tmp_path))


def test_esm2_token_budget_is_a_resource_not_a_credential(tmp_path):
    spec = JobSpec('esm2', ('--token-budget', '4096', '--api-key-env', 'PATHOGENWATCH_API_KEY'), tmp_path)
    assert '--token-budget' in spec.argv
