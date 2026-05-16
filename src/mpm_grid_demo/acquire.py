"""Benchmark acquisition helpers."""
from __future__ import annotations

from pathlib import Path
import subprocess

from .config import DemoConfig


def clone_andes_cases(config: DemoConfig) -> Path:
    """Clone the public ANDES cases repo if it is not already present."""
    dest = config.case_repo_dir
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists():
        return dest

    subprocess.run(["git", "clone", config.case_repo_url, str(dest)], check=True)
    return dest


def find_case_files(case_repo_dir: Path, keyword: str = "kundur") -> list[Path]:
    """Find plausible ANDES case files matching a keyword."""
    suffixes = {".xlsx", ".json", ".raw", ".dyr", ".m"}
    matches: list[Path] = []

    for path in case_repo_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in suffixes:
            continue
        if keyword.lower() in str(path).lower():
            matches.append(path)

    return sorted(matches)
