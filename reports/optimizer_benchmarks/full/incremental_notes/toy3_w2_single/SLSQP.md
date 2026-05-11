# Optimization Notes: toy3_w2_single / SLSQP

- Total runs so far: 4
- Feasible runs so far: 4
- Best feasible objective: 0.901532
- Best restart id: 2
- Best runtime (sec): 0.06101

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.905068 |     0.0959563 | Optimization terminated successfully |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.907043 |     0.0803377 | Optimization terminated successfully |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.901532 |     0.0610097 | Optimization terminated successfully |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.901532 |     0.0679822 | Optimization terminated successfully |