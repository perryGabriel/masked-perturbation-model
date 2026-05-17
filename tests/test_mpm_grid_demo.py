from __future__ import annotations

import importlib

import numpy as np
import control


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


def test_lft_builder_constructs_transfer_function_matrix():
    from mpmgame.lft import build_model_map

    A = np.array([[-1.0]])
    B_w = np.array([[2.0]])
    C_r = np.array([[3.0]])
    D_rw = np.array([[4.0]])

    M = build_model_map(A, B_w, C_r, D_rw)

    assert M.shape == (1, 1)
    assert isinstance(M[0, 0], control.TransferFunction)

    # M(s) = 3(s+1)^(-1)2 + 4 = 6/(s+1) + 4.
    value = control.evalfr(M[0, 0], 1j)
    expected = 6.0 / (1j + 1.0) + 4.0
    assert np.allclose(value, expected)


def test_bridge_constructs_lft_derived_transfer_function_matrix():
    from mpm_grid_demo.bridge import build_lft_from_linear_model
    from mpm_grid_demo.channels import select_random_channels
    from mpm_grid_demo.linearize import synthetic_two_area_like_model

    linear_model = synthetic_two_area_like_model()
    channels = select_random_channels(linear_model.A, linear_model.state_names, n_channels=2, random_seed=1)
    M = build_lft_from_linear_model(linear_model, channels)

    assert M.shape == (channels.C_r.shape[0], channels.B_w.shape[1])
    assert all(isinstance(M[i, j], control.TransferFunction) for i in range(M.shape[0]) for j in range(M.shape[1]))
