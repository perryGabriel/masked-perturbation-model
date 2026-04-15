import numpy as np
import mpmgame as mpm


def test_zero_sum_solver_on_reduced_game():
    payoff = np.array([[0.0, 1.0], [1.0, 0.0]])
    p, q, v = mpm.solve_zero_sum_game(payoff)
    np.testing.assert_allclose(p, [0.5, 0.5], atol=1e-7)
    np.testing.assert_allclose(q, [0.5, 0.5], atol=1e-7)
    assert abs(v - 0.5) < 1e-7
