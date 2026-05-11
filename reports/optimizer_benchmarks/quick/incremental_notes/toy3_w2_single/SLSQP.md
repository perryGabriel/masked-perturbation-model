# Optimization Notes: toy3_w2_single / SLSQP

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.897716
- Best restart id: 4
- Best runtime (sec): 0.1947

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.905068 |      0.215289 | Optimization terminated successfully |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.907044 |      0.111919 | Optimization terminated successfully |
|            2 |     33 | sparse_0         | True      | True       | False     |         0.90706  |      0.178114 | Optimization terminated successfully |
|            3 |     46 | gaussian_1       | True      | True       | False     |         0.897716 |      0.156132 | Optimization terminated successfully |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.897716 |      0.194739 | Optimization terminated successfully |