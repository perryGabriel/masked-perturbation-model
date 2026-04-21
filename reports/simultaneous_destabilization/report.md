# Simultaneous Destabilization Prototype Report

## Purpose
This report documents a finite-candidate simultaneous-destabilization workflow implemented in `mpmgame`.

The workflow analyzes:
- mask set `∇_1..∇_k`
- candidate attacks `Δ`
- per-mask destabilization success
- cross-success and overlap structure
- dominance relations and optional additive combinations

## Usage
Run the end-to-end demo:

```bash
python scripts/simultaneous_destabilization_demo.py
```

Run tests:

```bash
pytest tests/test_simultaneous_destabilization.py
```

Use from Python:

```python
from mpmgame import (
    AttackCandidate,
    MaskDefense,
    SimultaneousProblem,
    analyze_simultaneous_attacks,
)
```

## Assumptions
1. Finite attack candidate set (no continuous optimization over transfer functions).
2. Each mask and attack is represented on a shared finite support basis.
3. Default destabilization criterion is a weighted linear score exceeding a mask-specific threshold:
   `dot(mask * sensitivity, abs(attack)) >= threshold`.
4. Dominance is evaluated from cross-success row inclusion; ties can use `l1` norm.
5. Additive search is exploratory and combinatorial (default pairwise).

## Output schema
Output folder: `results/simultaneous_destabilization/`

### Machine-readable outputs
- `cross_success_matrix.csv`
  - rows: attack labels
  - columns: mask labels
  - values: `{0,1}` destabilization success
- `support_overlap_matrix.csv`
  - square mask-by-mask matrix
  - value: Jaccard overlap of support footprints
- `dominance_relations.json`
  - list of edges:
  - `[{"dominant": "Δi", "dominated": "Δj"}, ...]`
- `coverage_by_attack.json`
  - dictionary attack label -> list of mask labels destabilized
- `additive_attacks.json`
  - generated additive attack labels and vectors
- `run_summary.json`
  - aggregate counts and produced file paths

### Figures
- `cross_success_heatmap.png`
- `support_overlap_heatmap.png`
- `coverage_frequency.png`

## Key findings from the toy run
- Additive combinations expand the candidate set and can create broad-coverage attacks.
- Dominance edges are detectable even in small finite candidate pools.
- Support overlap helps explain when cross-mask success transfer is likely.
- The approach is deterministic for seeded toy problems.

## Reproduction commands
```bash
python scripts/simultaneous_destabilization_demo.py
python -m json.tool results/simultaneous_destabilization/run_summary.json
python -m json.tool results/simultaneous_destabilization/dominance_relations.json
```
