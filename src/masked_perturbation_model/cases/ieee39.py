"""IEEE 39-bus ANDES example case.

The module is safe to import without ANDES installed.  Loading/building the
actual ANDES system requires the optional ``andes`` extra.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

import numpy as np

from .andes.lft import AndesLFT, build_lft_from_state_matrix
from .andes.loaders import get_case_data_path, load_andes_case
from .andes.models import AndesCaseMetadata, AndesLinearizedSystem, AndesSystemModel

_METADATA_FILE = "ieee39_metadata.json"


def _load_metadata() -> AndesCaseMetadata:
    payload = json.loads(get_case_data_path(_METADATA_FILE).read_text(encoding="utf-8"))
    return AndesCaseMetadata(
        name=str(payload["name"]),
        description=str(payload["description"]),
        raw_case=str(payload["raw_case"]),
        dyr_case=payload.get("dyr_case"),
        expected_buses=payload.get("expected_buses"),
        expected_generators=payload.get("expected_generators"),
        source=str(payload.get("source", "ANDES built-in case library")),
        notes=tuple(str(note) for note in payload.get("notes", ())),
    )


def _get_state_names(system: Any, nstates: int) -> tuple[str, ...]:
    """Best-effort extraction of ANDES state names across ANDES versions."""

    dae = getattr(system, "dae", None)
    for attr in ("x_name", "x_name_output", "x_tex_name"):
        names = getattr(dae, attr, None)
        if names is not None:
            try:
                seq = tuple(str(item) for item in names)
            except TypeError:
                continue
            if len(seq) == nstates:
                return seq
    return tuple(f"x{i}" for i in range(nstates))


def _linearize(system: Any, metadata: AndesCaseMetadata) -> AndesLinearizedSystem:
    """Run ANDES' EIG small-signal routine and return state-space matrices."""

    eig = getattr(system, "EIG", None)
    if eig is None:
        raise RuntimeError("The loaded ANDES system does not expose an EIG routine")
    calc_as = getattr(eig, "calc_As", None)
    if callable(calc_as):
        A = calc_as(dense=True)
    else:
        run = getattr(eig, "run", None)
        if not callable(run):
            raise RuntimeError("The loaded ANDES system cannot calculate a state matrix")
        run()
        A = getattr(eig, "As", None)
    if A is None:
        A = getattr(eig, "As", None)
    if A is None:
        raise RuntimeError("ANDES did not return a state matrix from EIG")
    A_arr = np.asarray(A, dtype=float)
    if A_arr.ndim != 2 or A_arr.shape[0] != A_arr.shape[1]:
        raise RuntimeError(f"Expected a square ANDES state matrix, got shape {A_arr.shape}")
    n = int(A_arr.shape[0])
    names = _get_state_names(system, n)
    return AndesLinearizedSystem(
        A=A_arr,
        B=np.eye(n),
        C=np.eye(n),
        D=np.zeros((n, n)),
        state_names=names,
        input_names=tuple(f"u_{name}" for name in names),
        output_names=tuple(f"y_{name}" for name in names),
        metadata={"case": metadata.name, "source": metadata.source},
    )


@dataclass(frozen=True)
class IEEE39Model:
    """Loader and builder for the ANDES IEEE 39-bus example case."""

    metadata: AndesCaseMetadata

    @classmethod
    def from_package_data(cls) -> "IEEE39Model":
        """Create the IEEE39 case descriptor from packaged metadata."""

        return cls(metadata=_load_metadata())

    def load_andes_system(self, *, setup: bool = True, run_power_flow: bool = True, **kwargs: Any) -> Any:
        """Load the underlying ANDES system and optionally initialize it.

        This method requires the optional ANDES dependency.  Package imports and
        metadata inspection remain available without ANDES installed.
        """

        system = load_andes_case(self.metadata.raw_case, addfile=self.metadata.dyr_case, setup=False, **kwargs)
        if setup:
            setup_fn = getattr(system, "setup", None)
            if callable(setup_fn):
                setup_fn()
        if run_power_flow:
            pflow = getattr(system, "PFlow", None)
            run = getattr(pflow, "run", None)
            if callable(run):
                run()
        return system

    def build_system(
        self,
        *,
        setup: bool = True,
        run_power_flow: bool = True,
        linearize: bool = True,
        **kwargs: Any,
    ) -> AndesSystemModel:
        """Build the ANDES IEEE39 system and, by default, its state matrix."""

        system = self.load_andes_system(setup=setup, run_power_flow=run_power_flow, **kwargs)
        linearized = _linearize(system, self.metadata) if linearize else None
        return AndesSystemModel(metadata=self.metadata, andes_system=system, linearized=linearized)

    def build_lft(self, *, system_model: AndesSystemModel | None = None, **kwargs: Any) -> AndesLFT:
        """Build an LFT-style representation from the IEEE39 linearized model."""

        model = system_model if system_model is not None else self.build_system(**kwargs)
        if model.linearized is None:
            raise ValueError("system_model must include linearized state-space data")
        return build_lft_from_state_matrix(
            model.linearized.A,
            state_names=model.linearized.state_names,
            metadata={"case": self.metadata.name, "source": self.metadata.source},
        )

    def summary(self) -> dict[str, Any]:
        """Return static case metadata without importing ANDES."""

        return {
            "name": self.metadata.name,
            "description": self.metadata.description,
            "raw_case": self.metadata.raw_case,
            "dyr_case": self.metadata.dyr_case,
            "expected_buses": self.metadata.expected_buses,
            "expected_generators": self.metadata.expected_generators,
            "source": self.metadata.source,
        }


def load_ieee39_case() -> IEEE39Model:
    """Return the installable IEEE39 ANDES case descriptor."""

    return IEEE39Model.from_package_data()


def build_ieee39_lft(**kwargs: Any) -> AndesLFT:
    """Convenience wrapper for ``load_ieee39_case().build_lft(...)``."""

    return load_ieee39_case().build_lft(**kwargs)


__all__ = ["IEEE39Model", "build_ieee39_lft", "load_ieee39_case"]
