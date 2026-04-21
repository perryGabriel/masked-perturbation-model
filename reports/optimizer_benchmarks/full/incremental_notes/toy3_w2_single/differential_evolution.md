# Optimization Notes: toy3_w2_single / differential_evolution

- Total runs so far: 4
- Feasible runs so far: 2
- Best feasible objective: 0.910556
- Best restart id: 2
- Best runtime (sec): 1.825

## Failure modes observed

timeout: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                               |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |       inf        |       5       | timeout                               |
|            1 |     24 | gaussian_0       | False     | False      | True      |       inf        |       5       | timeout                               |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.910556 |       1.82526 | Optimization terminated successfully. |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.916985 |       1.54532 | Optimization terminated successfully. |