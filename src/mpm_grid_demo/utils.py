"""Small utility functions."""
from __future__ import annotations

from itertools import chain, combinations
from pathlib import Path
import subprocess
from typing import Iterable


def run(cmd: list[str], cwd: str | Path | None = None) -> None:
    print("$", " ".join(map(str, cmd)))
    subprocess.run(cmd, cwd=cwd, check=True)


def powerset(iterable: Iterable[int]):
    items = list(iterable)
    return chain.from_iterable(combinations(items, r) for r in range(len(items) + 1))


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
