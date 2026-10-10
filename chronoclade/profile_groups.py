"""Threshold sensitivity for local genetic groups, separate from time/place patterns."""
from collections import Counter
import numpy as np


def group_sensitivity(matrix, *, selected_threshold, candidates=None):
    """Summarise complete-linkage cuts, without treating a PCA separation as evidence.

    Thresholds are exploratory profile-distance fractions. No scheme-independent
    outbreak or clone threshold is inferred. The selected explicit threshold is
    included, and the same all-pairs diameter rule applies at every cut.
    """
    values = sorted(set(candidates or [selected_threshold / 2, selected_threshold,
                                      min(1, selected_threshold * 2)]))
    if any(not np.isfinite(value) or not 0 <= value <= 1 for value in values):
        raise ValueError('Group sensitivity thresholds must be finite fractions in [0,1]')
    result = []
    for cutoff in values:
        sizes = [len(group) for group in complete_linkage_groups(matrix, cutoff)]
        result.append(dict(distance_fraction=cutoff, group_count=len(sizes),
                           singletons=sizes.count(1), largest_group=max(sizes, default=0),
                           size_distribution=dict(sorted(Counter(sizes).items())),
                           selected=cutoff == selected_threshold))
    return dict(method='complete linkage; every within-group pair meets the cutoff',
                threshold_basis='explicit descriptive fraction of jointly called loci',
                thresholds=result,
                calibration_status='exploratory_not_biologically_validated',
                interpretation='Compare membership across thresholds and locus resamples; '
                'these groups do not establish rapid expansion, persistence, transmission or outbreaks.')


def complete_linkage_groups(matrix, threshold):
    """Deterministic complete linkage: all members satisfy the distance threshold."""
    groups = [(i,) for i in range(len(matrix))]
    linkage = matrix.copy()
    np.fill_diagonal(linkage, np.inf)
    while len(groups) > 1:
        position = int(np.argmin(linkage))
        i, j = np.unravel_index(position, linkage.shape)
        if linkage[i, j] > threshold:
            break
        if i > j:
            i, j = j, i
        groups[i] = tuple(sorted(groups[i] + groups[j]))
        groups.pop(j)
        linkage[i, :] = np.maximum(linkage[i, :], linkage[j, :])
        linkage[:, i] = linkage[i, :]
        linkage = np.delete(np.delete(linkage, j, axis=0), j, axis=1)
        np.fill_diagonal(linkage, np.inf)
    return sorted(groups)

