"""Import smoke tests for installable example cases."""


def test_masked_perturbation_model_cases_imports():
    from masked_perturbation_model.cases import ieee39
    from masked_perturbation_model.cases import load_ieee39_case
    from masked_perturbation_model.cases.ieee39 import IEEE39Model, build_ieee39_lft

    case = load_ieee39_case()
    assert isinstance(case, IEEE39Model)
    assert ieee39.load_ieee39_case().metadata.name == "ieee39"
    assert callable(build_ieee39_lft)


def test_mpmgame_cases_compatibility_imports():
    from mpmgame.cases import load_ieee39_case

    assert load_ieee39_case().summary()["expected_buses"] == 39


def test_package_import_does_not_require_andes():
    import masked_perturbation_model as mpm
    from masked_perturbation_model.cases import load_ieee39_case

    assert hasattr(mpm, "paper_example_data")
    assert load_ieee39_case().summary()["name"] == "ieee39"
