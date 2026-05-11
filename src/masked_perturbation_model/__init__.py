"""Compatibility package for the masked-perturbation model project.

The distribution is currently named ``mpmgame``.  This package provides the
underscore import path requested by downstream examples while preserving the
existing ``mpmgame`` public API.
"""

from mpmgame import *  # noqa: F401,F403

try:
    from mpmgame import __all__ as _mpmgame_all
except ImportError:  # pragma: no cover
    _mpmgame_all = []

__all__ = list(_mpmgame_all)
