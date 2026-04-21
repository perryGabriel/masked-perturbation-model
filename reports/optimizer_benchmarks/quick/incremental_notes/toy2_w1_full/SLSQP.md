# Optimization Notes: toy2_w1_full / SLSQP

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.702491
- Best restart id: 3
- Best runtime (sec): 0.002445

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.702491 |     0.0104081 | Optimization terminated successfully |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0018198 | Optimization terminated successfully |
|            2 |     33 | sparse_0         | True      | True       | False     |         0.702491 |     0.0026684 | Optimization terminated successfully |
|            3 |     46 | gaussian_1       | True      | True       | False     |         0.702491 |     0.0024447 | Optimization terminated successfully |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.702491 |     0.003355  | Optimization terminated successfully |