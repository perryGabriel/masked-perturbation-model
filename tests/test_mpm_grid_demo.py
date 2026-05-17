from __future__ import annotations

import importlib

import numpy as np


def test_mpm_grid_demo_package_can_be_imported():
    module = importlib.import_module("mpm_grid_demo")
    assert module.__version__


def test_synthetic_fallback_is_strictly_stable():
    from mpm_grid_demo.linearize import synthetic_two_area_like_model, verify_stability

    linear_model = synthetic_two_area_like_model()
    stable, eigvals = verify_stability(linear_model.A)

    assert stable
    assert eigvals.shape == (linear_model.A.shape[0],)


def test_random_channel_selection_is_dimensionally_consistent():
    from mpm_grid_demo.channels import select_random_channels
    from mpm_grid_demo.linearize import synthetic_two_area_like_model

    linear_model = synthetic_two_area_like_model()
    channels = select_random_channels(linear_model.A, linear_model.state_names, n_channels=2, random_seed=1)

    assert channels.B_w.shape == (linear_model.A.shape[0], 2)
    assert channels.C_r.shape == (2, linear_model.A.shape[0])
    assert channels.D_rw.shape == (2, 2)


def test_bridge_constructs_model_map_with_frequency_eval():
    from mpm_grid_demo.bridge import build_lft_from_linear_model
    from mpm_grid_demo.channels import select_random_channels
    from mpm_grid_demo.linearize import synthetic_two_area_like_model

    linear_model = synthetic_two_area_like_model()
    channels = select_random_channels(linear_model.A, linear_model.state_names, n_channels=2, random_seed=1)
    model_map = build_lft_from_linear_model(linear_model, channels)

    assert model_map.ninputs == channels.B_w.shape[1]
    assert model_map.noutputs == channels.C_r.shape[0]
    assert model_map.nstates == linear_model.A.shape[0]
    assert np.asarray(model_map.eval(1j)).shape == (channels.C_r.shape[0], channels.B_w.shape[1])
