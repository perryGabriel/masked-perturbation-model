import numpy as np
import mpmgame as mpm


def test_measured_respects_single_link_lower_bound():
    G = np.diag([0.5, 0.4, 0.3])
    alpha = np.array([[0.0, 0.1, 0.0], [0.05, 0.0, 0.1], [0.0, 0.05, 0.0]])
    Q = np.array([[0.0, 0.1, 0.0], [0.05, 0.0, 0.05], [0.0, 0.08, 0.0]])
    system = mpm.build_contract_system(G, alpha)
    pbar = mpm.access_matrix(system, Q, "w2")

    measured = mpm.vulnerability_single_link(system, Q, pbar).value
    lb = mpm.bound_lower_single_link(system, Q, pbar)
    assert measured + 1e-9 >= lb.value


def test_measured_not_above_upper_full_helper():
    G = np.diag([0.3, 0.2])
    alpha = np.array([[0.0, 0.05], [0.05, 0.0]])
    Q = np.zeros((2, 2))
    system = mpm.build_contract_system(G, alpha)
    pbar = mpm.access_matrix(system, Q, "w2")

    measured = mpm.vulnerability_full(system, Q, pbar).value
    ub = mpm.bound_upper_full(system, pbar)
    assert measured <= ub.value + 1e-9
