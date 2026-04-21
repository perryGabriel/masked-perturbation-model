# Optimization Notes: toy3_w2_single / SLSQP

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.897716
- Best restart id: 4
- Best runtime (sec): 0.05647

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.905068 |     0.0428425 | Optimization terminated successfully |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.907044 |     0.0226714 | Optimization terminated successfully |
|            2 |     33 | sparse_0         | True      | True       | False     |         0.90706  |     0.0432674 | Optimization terminated successfully |
|            3 |     46 | gaussian_1       | True      | True       | False     |         0.897716 |     0.0459505 | Optimization terminated successfully |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.897716 |     0.0564676 | Optimization terminated successfully |