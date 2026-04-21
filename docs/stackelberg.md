# Stackelberg masked game guide

This guide documents the finite-action Stackelberg utilities in `mpmgame.stackelberg` and the end-to-end demo in `scripts/stackelberg_demo.py`.

## What is implemented

- Defender commitment over budget-feasible masks.
- Attacker best response against defender pure/mixed commitments.
- Mixed and pure Stackelberg solvers with adversarial attacker tie-breaking.
- Strategy/export helpers for machine-readable output tables.

All utilities assume an attacker-utility matrix `U` with shape `(n_defenses, n_attacks)`, where lower values are better for the defender.

## CLI usage

Run the toy scenario and export result artifacts:

```bash
python scripts/stackelberg_demo.py
```

This creates:

- `results/stackelberg/` for CSV/JSON outcome files
- `reports/stackelberg/` for visual summaries

## Python usage

```python
import numpy as np
from mpmgame import (
    StackelbergMask,
    StackelbergAttack,
    StackelbergInstance,
    solve_stackelberg_mixed,
    solve_stackelberg_pure,
)

masks = [
    StackelbergMask("∇0", cost=0.0),
    StackelbergMask("∇1", cost=1.0),
    StackelbergMask("∇2", cost=2.0),
]
attacks = [StackelbergAttack("Δ1"), StackelbergAttack("Δ2")]
U = np.array([
    [0.9, 1.0],
    [0.6, 0.8],
    [0.5, 0.7],
], dtype=float)

instance = StackelbergInstance(payoff=U, masks=masks, attacks=attacks, budget=1.0)

mixed = solve_stackelberg_mixed(instance)
pure = solve_stackelberg_pure(instance)

print(mixed.defender_strategy, mixed.attacker_value)
print(pure.defender_strategy, pure.attacker_value)
```

## Output schema

### `results/stackelberg/*_summary.json`

- `mode`: Stackelberg mode (`"mixed"` or `"pure"`)
- `budget`: budget used for feasible mask filtering
- `feasible_indices`: defense row indices allowed by budget
- `attacker_value`: attacker best-response value at the commitment
- `attacker_best_response_indices`: attack indices in the best-response set
- `attacker_best_response_labels`: attack labels aligned to indices
- `defender_strategy`: full defender probability vector over all masks

### `results/stackelberg/*_policy_table.csv`

Columns:

- `mask_index`
- `mask_label`
- `mask_cost`
- `feasible`
- `probability`

### Additional demo outputs

- `toy_payoff_matrix.csv`: demo payoff matrix with defense metadata.
- `mixed_expected_attack_utilities.csv`: expected attacker utility per attack under mixed commitment.
- `reports/stackelberg/response_map.png`: pure-response map by defender action.
- `reports/stackelberg/value_comparison.png`: mixed vs pure worst-case attacker value.

## Caveats

- This module solves finite-action commitment only; it does not synthesize new masks.
- The default objective is adversarial tie-breaking (attacker chooses the worst defender outcome among ties).
- For zero-sum finite games, Stackelberg mixed value can match minimax value; interpret differences with care.
- Budget handling is action filtering (`cost <= budget`), not a knapsack over simultaneously deployed masks.
