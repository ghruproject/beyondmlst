"""The optional exploration route shares profile-only context acquisition."""
from pathlib import Path
from types import SimpleNamespace
from typer.testing import CliRunner
from chronoclade.cli import app


def test_genome_context_handoff_uses_fetched_dataset(monkeypatch, tmp_path):
    calls = []
    original = tmp_path / 'queries.json'
    retrieved = tmp_path / 'with-context/dataset.json'
    def discover(source, output, **kwargs):
        calls.append((source, output, kwargs))
        return SimpleNamespace(dataset_manifest=retrieved)
    def analyse(source, output, **kwargs):
        assert source == retrieved
        assert kwargs['catalogue_manifest'] == Path('alleles.json')
        assert kwargs['lin_level'] == 6
        return SimpleNamespace(manifest_path=output/'esm2.json', report_path=output/'report.html')
    monkeypatch.setattr('chronoclade.profile_context_provider.discover_profile_context', discover)
    monkeypatch.setattr('chronoclade.esm2.workflow.run_esm2', analyse)
    result = CliRunner().invoke(app, ['esm2', str(original), '--fetch-context', '--lin-level','6',
                                    '--allele-catalogue','alleles.json', '--out',str(tmp_path/'esm2')])
    assert result.exit_code == 0, result.output
    assert calls[0][0] == original
    assert calls[0][1] == tmp_path/'esm2_context'
    assert calls[0][2]['lin_level'] == 6


def test_profile_context_requires_genome_input():
    result = CliRunner().invoke(app, ['esm2','proteins.faa','--fetch-context'])
    assert result.exit_code == 2
    assert 'prepared-genome exploration' in result.output
