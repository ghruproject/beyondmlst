import numpy as np
from chronoclade.profile_groups import group_sensitivity


def test_groups_use_diameter_rule_and_explicit_threshold():
    matrix = np.array([[0, .01, .05], [.01, 0, .04], [.05, .04, 0]])
    result = group_sensitivity(matrix, selected_threshold=.02)
    assert result['calibration_status'] == 'exploratory_not_biologically_validated'
    selected = next(row for row in result['thresholds'] if row['selected'])
    assert selected['group_count'] == 2 and selected['largest_group'] == 2
