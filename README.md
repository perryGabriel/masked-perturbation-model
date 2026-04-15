# mpmgame: masked-perturbation security games for feedback systems

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
