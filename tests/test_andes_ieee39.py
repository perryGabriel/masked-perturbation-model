"""Tests for the ANDES IEEE39 example-system API."""

import numpy as np
import pytest

from masked_perturbation_model.cases import load_ieee39_case
from masked_perturbation_model.cases.andes.lft import build_lft_from_state_matrix
from masked_perturbation_model.cases.andes.loaders import require_andes


def test_ieee39_metadata_available_without_andes():
    case = load_ieee39_case()
    summary = case.summary()

    assert summary["name"] == "ieee39"
    assert summary["expected_buses"] == 39
    assert summary["expected_generators"] == 10
    assert summary["raw_case"].endswith("ieee39.raw")
    assert summary["dyr_case"].endswith("ieee39.dyr")


def test_optional_andes_dependency_error_is_helpful():
    try:
        require_andes()
    except ImportError as exc:
        assert "pip install" in str(exc)
        assert "mpmgame[andes]" in str(exc)
    else:
        pytest.skip("ANDES is installed; missing-dependency error path is not applicable")


def test_lft_container_dimensions_from_state_matrix():
    A = np.array([[-1.0, 0.25], [0.0, -2.0]])
    lft = build_lft_from_state_matrix(A, state_names=("delta", "omega"), metadata={"case": "toy"})

    assert lft.nstates == 2
    assert lft.ninputs == 2
    assert lft.noutputs == 2
    assert lft.A.shape == (2, 2)
    assert lft.B.shape == (2, 2)
    assert lft.C.shape == (2, 2)
    assert lft.D.shape == (2, 2)
    assert lft.mask_shape == (2, 2)
    assert lft.state_names == ("delta", "omega")


def test_ieee39_andes_build_smoke_when_dependency_available():
    pytest.importorskip("andes")

    case = load_ieee39_case()
    system = case.build_system()
    assert system.metadata.name == "ieee39"
    assert system.linearized is not None
    assert system.linearized.A.shape[0] == system.linearized.A.shape[1]
    assert system.linearized.nstates > 0

    lft = case.build_lft(system_model=system)
    assert lft.A.shape == system.linearized.A.shape
    assert lft.B.shape[0] == lft.nstates
    assert lft.C.shape[1] == lft.nstates
