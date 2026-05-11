"""Loaders for optional ANDES example-system dependencies and package data."""

from __future__ import annotations

from importlib import import_module
from importlib.resources import files
from pathlib import Path
from types import ModuleType
from typing import Any

DISTRIBUTION_NAME = "mpmgame"
ANDES_EXTRA = f'{DISTRIBUTION_NAME}[andes]'
DATA_PACKAGE = "masked_perturbation_model.cases.andes.data"


def require_andes() -> ModuleType:
    """Return the ANDES module or raise a helpful optional-dependency error."""

    try:
        return import_module("andes")
    except ImportError as exc:  # pragma: no cover - exact message tested indirectly
        raise ImportError(
            "ANDES support is optional. Install it with "
            f'`pip install "{ANDES_EXTRA}"` before loading ANDES-backed cases.'
        ) from exc


def get_case_data_path(filename: str) -> Path:
    """Return the path for a packaged ANDES case metadata/data file."""

    return Path(str(files(DATA_PACKAGE).joinpath(filename)))


def resolve_andes_case(case_name: str) -> str:
    """Resolve an ANDES built-in case name when ANDES is installed.

    The function first asks ``andes.get_case`` for the packaged case path and
    falls back to the provided name to support user installations that expose
    cases directly on the ANDES search path.
    """

    andes = require_andes()
    get_case = getattr(andes, "get_case", None)
    if callable(get_case):
        try:
            return str(get_case(case_name))
        except Exception:
            return case_name
    return case_name


def load_andes_case(case_name: str, *, addfile: str | None = None, setup: bool = False, **kwargs: Any) -> Any:
    """Load an ANDES case using the optional ANDES dependency.

    Parameters mirror ``andes.load`` for the subset used by the example cases.
    ``case_name`` and ``addfile`` may be ANDES built-in case names such as
    ``ieee39/ieee39.raw`` and ``ieee39/ieee39.dyr``.
    """

    andes = require_andes()
    load = getattr(andes, "load")
    resolved_case = resolve_andes_case(case_name)
    resolved_addfile = resolve_andes_case(addfile) if addfile else None
    if resolved_addfile is None:
        return load(resolved_case, setup=setup, **kwargs)
    return load(resolved_case, addfile=resolved_addfile, setup=setup, **kwargs)
