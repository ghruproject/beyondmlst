"""Binary profile distances preserve selection ties and portable evidence."""
import shutil
import numpy as np
import pytest

from chronoclade.cgmlst.scale_distance import write_distances
from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import load_selection_ensemble, load_selection_manifest, write_selection_ensemble


def test_binary_selection_roundtrip_and_tampering(tmp_path):
    rows = [dict(sample_id=ident, origin='local' if ident == 'q' else 'context',
                 role='input' if ident == 'q' else 'context', country='Greece',
                 species='Klebsiella pneumoniae', cgmlst_scheme='test', cgmlst_scheme_version='1',
                 cgmlst_loci=['a','b','c','d'], cgmlst_profile=dict(zip(['a','b','c','d'], calls)))
            for ident, calls in [('q',[1,1,1,1]), ('x',[2,1,1,1]), ('y',[2,1,1,0]), ('far',[2,2,2,2])]]
    evidence, _, _ = write_distances(rows, tmp_path / 'matrix', min_overlap=.7)
    prepared = write_dataset(PreparedDataset(samples=samples_from_records(rows)), tmp_path / 'prepared')
    manifest = write_selection_ensemble(rows[:1], rows[1:], {'distance_evidence':evidence},
                                        dataset_manifest=prepared, output=tmp_path / 'runs',
                                        size=0, nearest_per_query=1, replicates=2)
    ensemble = load_selection_ensemble(manifest)
    selection = load_selection_manifest(manifest.parent / ensemble['selections'][0]['path'])
    assert set(selection['selected_sample_ids']) == {'q','x','y'}
    assert selection['budget_overrun'] == 2
    moved = tmp_path / 'moved'
    moved.mkdir()
    for name in ('prepared','runs'):
        shutil.move(str(tmp_path / name), moved / name)
    shutil.rmtree(tmp_path / 'matrix')
    restored = load_selection_ensemble(moved / 'runs/ensemble.json')
    assert restored == ensemble
    matrix = np.load(moved / 'runs/distance.npy', mmap_mode='r+')
    matrix[0,1] = 0
    matrix.flush()
    with pytest.raises((ValueError, WorkflowError), match='checksum'):
        load_selection_ensemble(moved / 'runs/ensemble.json')
