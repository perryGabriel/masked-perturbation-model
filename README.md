# mpmgame: masked-perturbation security games for feedback systems

A YouTube video of the defense of this Master's Thesis work is available [here](https://youtu.be/jjW-QZsh_lg).

`mpmgame` is a focused research package that reproduces the masked-perturbation model (MPM) finite-action game formulation at the transfer-function level.

It models:
- feedback map `M`
- attack operator `Δ`
- binary defense mask `∇`

with core equations:

- `r = M w + d1`
- `w = (Δ ∘ ∇) r + d2`

A `0` in `∇` blocks a link; a `1` leaves it available to the attacker.

## Scope

This package is for finite-action masked-perturbation games and paper-style worked examples.

It is **not** a full robust-control toolbox and does **not** implement structured singular value (μ) analysis.

## Installation

```bash
pip install -e .
```

## Quick start

```python
import mpmgame as mpm

data = mpm.paper_example_data()
result = mpm.run_paper_example()

print(result.success_sets.attack_success)
print(result.reduced.payoff)
print(result.attacker_mix, result.defender_mix, result.value)
```

## Core features

- Rank-1 defense masks and defense costs
- Admissible-defense enumeration under budget
- Mask algebra (`union`, `intersection`, `complement`, subset relation)
- Destabilization checks from transfer-function feedback interconnections
- Success sets and dominated-strategy elimination
- Zero-sum mixed equilibrium via linear programming (`scipy.optimize.linprog`)
- Reproduction of the paper’s 2x2 reduced-game equilibrium


## Python package module guide

The Python implementation lives under `src/mpmgame/`:

- `core.py`: Main user-facing API for masks, costs, masking operations, stability utilities, success sets, dominance, and payoff construction.
- `tf_tools.py`: Transfer-function matrix helpers and closed-loop pole/stability checks.
- `game.py`: Zero-sum finite-game LP solver and expected utility.
- `examples.py`: Paper-example constructors (`paper_example_data`, `run_paper_example`) and coordinated attack helper.
- `plotting.py`: Optional plotting helpers (payoff heatmap and success bipartite graph).

For quick interactive use, import directly from the package root (`import mpmgame as mpm`) since these APIs are re-exported in `mpmgame.__init__`.

### DSF naming and migration note

- DSF-specific exports are now available under explicit names:
  - `DynamicQ_dsf`
  - `serialize_dynamic_q_for_row_dsf(...)`
  - `deserialize_dynamic_q_from_row_dsf(...)`
  - `benchmark_problem_registry_dsf(...)` (opt-in DSF toy registry)
- The default `benchmark_problem_registry(...)` remains static (`ParameterizationSpec(kind="static_hollow")`) for reproducibility with existing notebooks/scripts.
- Internal opt-in constructors now include DSF variants (for example `_toy_problem_dsf(...)`) rather than silently changing static defaults.
- **Migration note:** “No behavior change unless `ParameterizationSpec.kind='dsf_poly'`”.

## Repository layout

```text
src/mpmgame/      Python package
tests/            Pytest suite
notebooks/        Replication notebook
matlab/           Preserved MATLAB material
docs/             Notes
```

## Run tests

```bash
pytest
```

## Notebook

See `notebooks/paper_example.ipynb` for an end-to-end replication of the paper example.


## Chapter-2 FMP extension (numerical demos)

This repo now includes chapter-2 feedback-vulnerability utilities for small-system numerical exploration under explicit threat/access assumptions.

Added modules:
- `mpmgame.fmp`: contract-system setup (`Gamma`, `P`, attack-map) and measured vulnerabilities.
- `mpmgame.bounds`: theorem-inspired bound helpers (with assumptions documented in docstrings).
- `mpmgame.experiments_ch2`: family sweeps (`vary_q`, `vary_alpha`, `vary_access`) and report export.
- `mpmgame.random_systems`: conservative random-system filtering and bound-ratio data.
- `mpmgame.design_projection`: projection-inspired finite-dimensional surrogate iterations.
- `mpmgame.design_lp`: LP-inspired finite-dimensional relaxation via `scipy.optimize.linprog`.
- `mpmgame.plotting_ch2`: plotting helpers for vulnerability-vs-parameter and design diagnostics.
- `mpmgame.realization_report`: helpers to turn optimizer CSV outputs into per-result markdown notes and inspect the recovered realization (`Q`, `P=(I-Q)G`) in numeric/SymPy-friendly form.

### Exact vs approximate

- **Exact for static matrices:** `Gamma`, `P=(I-Q)G`, measured full/single-link vulnerability, and matrix-norm calculations.
- **Approximate / surrogate:** projection and LP design methods are finite-dimensional static relaxations inspired by the chapter notes.
- **Not a proof engine:** outputs are numerical sanity checks on toy systems only.

### Chapter-2 notebooks

- `notebooks/chapter2_bounds_demo.ipynb`
- `notebooks/chapter2_design_demo.ipynb`

### Chapter-2 docs

- `docs/chapter2.md` describes assumptions and interpretation.

### Nonlinear optimizer reporting workflow

To keep per-result interpretation notes while runs are still in progress, pass `incremental_notes_dir` to `run_benchmark_suite(...)`. This creates one markdown file per `(problem_id, algorithm)` and updates it after each restart.

After a run, you can emit a realization-inspection markdown document from the raw CSV:

```python
from mpmgame.realization_report import write_realization_markdown_from_csv

write_realization_markdown_from_csv(
    "results/optimizer_benchmarks/full/benchmark_raw_results.csv",
    "reports/optimizer_benchmarks/full/realizations/toy3_w3_full_best.md",
    problem_id="toy3_w3_full",
)
```
