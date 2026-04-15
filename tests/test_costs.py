import numpy as np
import mpmgame as mpm


def test_defense_cost():
    m = np.array([[1, 0], [0, 0]])
    c = mpm.defense_cost(mask=m, c_w=[2, 1], c_r=[3, 1])
    assert c == 2.0


def test_admissible_defenses_example_budget():
    defs = mpm.admissible_defenses(2, 2, c_w=[2, 1], c_r=[3, 1], budget=2, include_empty=False)
    expected = [
        np.array([[1, 0], [1, 0]]),
        np.array([[1, 1], [0, 0]]),
        np.array([[0, 0], [1, 1]]),
        np.array([[1, 0], [0, 0]]),
    ]
    got = {tuple(m.flatten()) for m in defs}
    want = {tuple(m.flatten()) for m in expected}
    assert got == want
