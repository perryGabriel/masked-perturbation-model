import mpmgame as mpm


def test_success_sets_match_paper_example():
    data = mpm.paper_example_data()
    s = mpm.compute_success_sets(data.M, data.attacks, data.defenses)

    assert s.attack_success["Δ1"] == {"∇1", "∇2", "∇4"}
    assert s.attack_success["Δ2"] == {"∇3"}
    assert s.attack_success["Δ3"] == {"∇2"}

    assert s.defense_success["∇1"] == {"Δ2", "Δ3"}
    assert s.defense_success["∇2"] == {"Δ2"}
    assert s.defense_success["∇3"] == {"Δ1", "Δ3"}
    assert s.defense_success["∇4"] == {"Δ2", "Δ3"}
