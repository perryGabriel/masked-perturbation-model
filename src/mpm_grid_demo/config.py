"""Configuration for the ANDES-backed MPM verification demo."""
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DemoConfig:
    project_root: Path = Path(__file__).resolve().parents[2]
    case_repo_url: str = "https://github.com/CURENT2/andes_cases.git"
    case_repo_dir: Path = project_root / "data" / "raw" / "andes_cases"
    preferred_case_keyword: str = "ieee39"
    defense_budget: float = 2.0
    top_k_single_link_attacks: int = 5
    frequency_min: float = 1e-3
    frequency_max: float = 1e3
    frequency_points: int = 2000
    attack_gain_safety_factor: float = 1.05
    stability_tol: float = 1e-8
    use_synthetic_fallback: bool = True
    random_seed: int = 42
