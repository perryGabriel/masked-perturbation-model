import mpmgame as mpm


def test_dominance_and_reduced_game():
    data = mpm.paper_example_data()
    ss = mpm.compute_success_sets(data.M, data.attacks, data.defenses)
    datk = mpm.dominated_attacks(ss.attack_success)
    ddef = mpm.dominated_defenses(ss.defense_success)

    assert "Δ3" in datk
    assert "∇1" in ddef
    assert "∇2" in ddef

    reduced = mpm.eliminate_dominated_strategies(data.M, data.attacks, data.defenses)
    assert [a.label for a in reduced.attacks] == ["Δ1", "Δ2"]
    assert [d.label for d in reduced.defenses] == ["∇3", "∇4"]
