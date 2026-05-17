"""ANDES loading and initialization helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def load_andes_system(case_path: str | Path, setup: bool = False) -> Any:
    """Load an ANDES case.

    ANDES is imported lazily so that this demo package can be imported in
    lightweight environments where the optional dependency is not installed.
    """
    import andes

    return andes.load(str(case_path), setup=setup)


def setup_power_flow_and_dynamics(system: Any) -> Any:
    """Set up the system, solve power flow, and initialize dynamic machinery."""
    system.setup()

    if hasattr(system, "PFlow"):
        system.PFlow.run()

    if hasattr(system, "TDS") and hasattr(system.TDS, "init"):
        system.TDS.init()

    return system
