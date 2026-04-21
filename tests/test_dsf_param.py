from __future__ import annotations

import numpy as np

from mpmgame.dsf_param import DynamicQ, deserialize_dynamic_q_from_row, serialize_dynamic_q_for_row


def test_dynamic_q_enforces_hollow_and_mask() -> None:
    num = np.zeros((3, 3, 2), dtype=float)
    num[0, 0, 0] = 1.0
    num[0, 2, 1] = 2.0
    mask = np.ones((3, 3), dtype=bool)
    mask[0, 2] = False
    q = DynamicQ.from_tensors(num, mask=mask, zero_diagonal=True, enforce=True)

    assert np.allclose(q.num_coeffs[np.arange(3), np.arange(3), :], 0.0)
    assert np.allclose(q.num_coeffs[0, 2, :], 0.0)


def test_dynamic_q_to_tf_matrix_builds_object_matrix() -> None:
    num = np.zeros((2, 2, 2), dtype=float)
    num[0, 1, 0] = 0.5
    num[0, 1, 1] = 0.25
    den = np.array([0.2], dtype=float)
    q = DynamicQ.from_tensors(num, den_coeffs=den)

    tfm = q.to_tf_matrix()
    assert tfm.shape == (2, 2)
    # off-diagonal should retain a first-order numerator and denominator.
    assert len(tfm[0, 1].num[0][0]) == 2
    assert len(tfm[0, 1].den[0][0]) == 2


def test_dynamic_q_payload_roundtrip() -> None:
    num = np.random.default_rng(0).normal(size=(2, 2, 3))
    den = np.random.default_rng(1).normal(size=(2, 2, 2))
    mask = np.array([[1, 1], [0, 1]], dtype=bool)
    q = DynamicQ.from_tensors(num, den_coeffs=den, mask=mask, metadata={"problem_id": "p1"})

    payload = serialize_dynamic_q_for_row(q)
    loaded = deserialize_dynamic_q_from_row(payload)
    assert loaded.shape == (2, 2)
    assert loaded.den_degree == 2
    assert loaded.metadata["problem_id"] == "p1"
