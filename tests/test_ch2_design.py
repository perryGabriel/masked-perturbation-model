import numpy as np
import mpmgame as mpm


def test_projection_reduces_surrogate_once():
    G = np.diag([0.7, 0.5, 0.4])
    alpha = np.array([[0.0, 0.2, 0.0], [0.1, 0.0, 0.1], [0.0, 0.1, 0.0]])
    system = mpm.build_contract_system(G, alpha)
    res = mpm.projection_design(system, max_iter=3, threat_model="single_link", access_model="w2")
    assert len(res.iterations) >= 1
    assert res.iterations[-1].surrogate_obj <= res.iterations[0].surrogate_obj + 1e-8


def test_lp_relaxation_feasible_on_toy():
    G = np.diag([0.6, 0.5, 0.4])
    alpha = np.array([[0.0, 0.1, 0.0], [0.1, 0.0, 0.1], [0.0, 0.05, 0.0]])
    system = mpm.build_contract_system(G, alpha)
    out = mpm.lp_relaxation_design(system, access_model="w2", threat_model="single_link", g_hat=0.2)
    assert out.success
