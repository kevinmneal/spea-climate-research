"""Small helpers for nested geographic selection; never accept outer scores."""
import numpy as np


def nested_roles(outer_train, inner_train, inner_test):
    """Intersect independently guarded roles with the outer training cohort."""
    arrays = [np.asarray(x, bool) for x in (outer_train, inner_train, inner_test)]
    if any(x.ndim != 1 or x.shape != arrays[0].shape for x in arrays):
        raise ValueError('Aligned one-dimensional role masks required')
    outer, train, test = arrays
    if np.any(train & test):
        raise ValueError('Inner train and validation roles overlap')
    return outer & train, outer & test


def select_candidate(scores, candidate_ids, minimum_folds=2, tolerance=1e-10):
    """Equal-fold mean; explicit candidate order breaks near-ties.

    Rows are candidates, columns are predeclared eligible inner folds. A failed
    candidate is not permitted to benefit from selective omission of folds.
    Empty/insufficient folds retain the first (fixed baseline) candidate.
    """
    values = np.asarray(scores, float)
    if values.ndim != 2 or values.shape[0] != len(candidate_ids) or not candidate_ids:
        raise ValueError('Candidate-by-inner-fold matrix required')
    if len(set(candidate_ids)) != len(candidate_ids) or not np.isfinite(values).all():
        raise ValueError('Unique candidates and complete finite inner scores required')
    if values.shape[1] < minimum_folds:
        return candidate_ids[0], 'insufficient_inner_folds_keep_baseline'
    means = values.mean(axis=1)
    chosen = np.flatnonzero(means >= means.max() - tolerance)[0]
    return candidate_ids[chosen], 'maximum_equal_inner_fold_block_score'
