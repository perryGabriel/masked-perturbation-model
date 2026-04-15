import numpy as np
import mpmgame as mpm


def test_gamma_and_realization_static():
    G = np.diag([0.6, 0.4])
    alpha = np.array([[0.0, 0.2], [0.1, 0.0]])
    Q = np.array([[0.0, 0.3], [0.2, 0.0]])
    system = mpm.build_contract_system(G, alpha)

    P = mpm.make_realization(G, Q)
    Gamma = mpm.compute_gamma(system.G, system.alpha)

    np.testing.assert_allclose(P, (np.eye(2) - Q) @ G)
    np.testing.assert_allclose(Gamma @ (np.eye(2) - G @ alpha), np.eye(2), atol=1e-8)


def test_access_models_w1_w2_w3():
    G = np.diag([0.8, 0.5])
    alpha = np.zeros((2, 2))
    Q = np.array([[0.0, 0.2], [0.1, 0.0]])
    system = mpm.build_contract_system(G, alpha)

    w1 = mpm.access_matrix(system, Q, "w1")
    w2 = mpm.access_matrix(system, Q, "w2")
    w3 = mpm.access_matrix(system, Q, "w3")

    np.testing.assert_allclose(w1, (np.eye(2) - Q) @ G)
    np.testing.assert_allclose(w2, np.eye(2))
    np.testing.assert_allclose(w3, Q)


def test_vulnerabilities_differ_on_unequal_structure():
    G = np.diag([0.7, 0.2])
    alpha = np.array([[0.0, 0.1], [0.0, 0.0]])
    Q = np.array([[0.0, 0.35], [0.02, 0.0]])
    system = mpm.build_contract_system(G, alpha)
    Pbar = mpm.access_matrix(system, Q, "w2")

    vf = mpm.vulnerability_full(system, Q, Pbar).value
    vs = mpm.vulnerability_single_link(system, Q, Pbar).value
    assert vf >= vs
    assert abs(vf - vs) > 1e-6
