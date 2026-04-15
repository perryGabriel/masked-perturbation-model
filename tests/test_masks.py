import numpy as np
import mpmgame as mpm


def test_rank1_mask():
    eta_w = [1, 0]
    eta_r = [1, 0]
    m = mpm.rank1_mask(eta_w, eta_r)
    np.testing.assert_array_equal(m, np.array([[1, 0], [0, 0]]))


def test_mask_set_operations():
    m1 = np.array([[1, 0], [1, 0]])
    m2 = np.array([[1, 1], [0, 0]])
    np.testing.assert_array_equal(mpm.union_masks(m1, m2), np.array([[1, 0], [0, 0]]))
    np.testing.assert_array_equal(mpm.complement_mask(m1), np.array([[0, 1], [0, 1]]))
    inter = mpm.intersect_masks(m1, m2)
    np.testing.assert_array_equal(inter, np.array([[1, 1], [1, 1]]))
    assert mpm.is_subset_mask(np.array([[1, 0], [1, 0]]), np.array([[1, 0], [0, 0]]))


def test_mask_footprint_counts():
    m = np.array([[1, 0], [0, 0]])
    fp = mpm.mask_footprint(m)
    assert fp["num_defended"] == 3
    assert fp["num_uncovered"] == 1
