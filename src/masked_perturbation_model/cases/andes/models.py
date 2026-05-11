"""Reusable data models for ANDES-backed example systems."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class AndesCaseMetadata:
    """Static metadata describing an ANDES example case."""

    name: str
    description: str
    raw_case: str
    dyr_case: str | None = None
    expected_buses: int | None = None
    expected_generators: int | None = None
    source: str = "ANDES built-in case library"
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class AndesLinearizedSystem:
    """Linearized state-space data extracted from an ANDES system.

    ``A`` is the small-signal state matrix returned by ANDES' EIG routine.  The
    default ``B``, ``C``, and ``D`` expose an identity state-input/state-output
    system so downstream MPM examples can inspect dimensions and channel names
    without inventing new control-design algorithms.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    state_names: tuple[str, ...] = ()
    input_names: tuple[str, ...] = ()
    output_names: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def nstates(self) -> int:
        return int(self.A.shape[0])

    @property
    def ninputs(self) -> int:
        return int(self.B.shape[1])

    @property
    def noutputs(self) -> int:
        return int(self.C.shape[0])


@dataclass
class AndesSystemModel:
    """Container returned by case builders for an initialized ANDES system."""

    metadata: AndesCaseMetadata
    andes_system: Any
    linearized: AndesLinearizedSystem | None = None

    @property
    def nstates(self) -> int | None:
        return None if self.linearized is None else self.linearized.nstates

    def summary(self) -> dict[str, Any]:
        """Return a lightweight summary safe for examples and smoke tests."""

        return {
            "name": self.metadata.name,
            "description": self.metadata.description,
            "expected_buses": self.metadata.expected_buses,
            "expected_generators": self.metadata.expected_generators,
            "nstates": self.nstates,
            "has_andes_system": self.andes_system is not None,
        }
