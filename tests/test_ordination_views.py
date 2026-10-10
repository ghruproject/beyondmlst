from pathlib import Path

import numpy as np
import pytest

from chronoclade.report_components.ordination import ordination_viewer, write_ordination_views


def test_shared_coordinates_metadata_views_keep_missing_dates(tmp_path):
    records = [
        dict(sample_id='input', role='input', country='Greece', collection_date='2020'),
        dict(sample_id='context', role='context', country='Italy', date_start='2021-01-01',
             date_end='2022-12-31', date_precision='interval'),
        dict(sample_id='missing', role='context'),
    ]
    paths = write_ordination_views(tmp_path, 'cohort', np.array([[0, 0], [1, 2], [3, 2]]), records)
    assert len(paths) == 4
    assert all(Path(path).exists() for path in paths.values())
    html = Path(paths['pcoa_views_html']).read_text()
    assert 'Colour by' in html and 'Country' in html and 'Date' in html
    date = Path(paths['pcoa_date_figure']).read_text()
    assert 'Date unknown' in date and 'Input genomes' in date and 'Public comparisons' in date
    assert 'iframe' in ordination_viewer([paths], tmp_path)
    assert ordination_viewer([dict(pcoa_views_html='../escape.html')], tmp_path) == ''


def test_invalid_coordinates_and_names_rejected(tmp_path):
    with pytest.raises(ValueError, match='two finite'):
        write_ordination_views(tmp_path, 'safe', [[np.nan, 0]], [dict(sample_id='a')])
    with pytest.raises(ValueError, match='safe artifact'):
        write_ordination_views(tmp_path, '../bad', [[0, 0]], [dict(sample_id='a')])
