import numpy as np
import mpmgame as mpm


def test_full_paper_example_outputs():
    result = mpm.run_paper_example()
    np.testing.assert_array_equal(result.reduced.payoff, np.array([[0.0, 1.0], [1.0, 0.0]]))
    np.testing.assert_allclose(result.attacker_mix, [0.5, 0.5], atol=1e-7)
    np.testing.assert_allclose(result.defender_mix, [0.5, 0.5], atol=1e-7)
    assert abs(result.value - 0.5) < 1e-7


from mpmgame.examples import coordinated_attack_defeats_all


def test_coordinated_attack_defeats_all():
    assert coordinated_attack_defeats_all()
