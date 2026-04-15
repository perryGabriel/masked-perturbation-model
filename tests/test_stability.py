import mpmgame as mpm


def test_masked_delta_and_masked_model_view_consistent():
    data = mpm.paper_example_data()
    for atk in data.attacks:
        for dfn in data.defenses:
            a = mpm.is_destabilizing(data.M, atk.delta, dfn.mask, method="mask_delta")
            b = mpm.is_destabilizing(data.M, atk.delta, dfn.mask, method="mask_model_map")
            assert a == b
