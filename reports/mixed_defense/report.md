# Mixed-defense action-set benchmark report

## Overview
This report documents the mixed-defense action-set prototype implemented in:

- `src/mpmgame/mixed_defense.py`
- `scripts/mixed_defense_benchmark.py`

The implementation targets roadmap section **2.2** by constructing capped defense subsets `|ζ| ≤ z` and evaluating reduced games over the selected defenses.

## Usage
Run the benchmark sweep:

```bash
python scripts/mixed_defense_benchmark.py
```

This runs parameter sweeps over:

- cardinality cap `z ∈ {1,2,3,4}`
- objective weights `(α_union, α_intersection, α_cardinality)`
- selection method (`greedy`, `random`)

## Expected outputs
### Result artifacts (`results/mixed_defense/`)
- `sweep_results.csv`: one row per `(z, weights, method)` trial with selection and reduced-game metrics.
- `summary.json`: full-game reference value and best rows per method.

### Figure artifacts (`reports/mixed_defense/`)
- `value_vs_z.png`: reduced-game value versus action-set size cap.
- `coverage_vs_z.png`: attack coverage of selected defenses versus `z`.
- `objective_vs_value_scatter.png`: relationship between surrogate objective and reduced-game value.

## Approximation-quality metrics
The benchmark records quality with:

- `value_gap = reduced_value - full_value` (distance from full defense set value)
- `union_coverage` (number of attacks blocked by at least one selected defense)
- `selected_size` (must satisfy `selected_size ≤ z`)

These metrics enable practical validation of reduced-set quality without claiming formal guarantees for all weight settings.

## Code and artifact paths
- API implementation: `src/mpmgame/mixed_defense.py`
- Reduced-game utilities: `src/mpmgame/core.py`
- Sweep runner: `scripts/mixed_defense_benchmark.py`
- Results directory: `results/mixed_defense/`
- Report/plots directory: `reports/mixed_defense/`
