# Optimization Notes: toy3_w3_full / SLSQP

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 2.51755e-06
- Best restart id: 1
- Best runtime (sec): 0.3926

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |      3.00447e-05 |     0.0270083 | Optimization terminated successfully |
|            1 |     20 | gaussian_0       | True      | True       | False     |      2.51755e-06 |     0.392587  | Optimization terminated successfully |
|            2 |     33 | sparse_0         | True      | True       | False     |      3.13308e-06 |     0.309514  | Optimization terminated successfully |
|            3 |     46 | gaussian_1       | True      | True       | False     |      2.55026e-05 |     0.291787  | Optimization terminated successfully |
|            4 |     59 | sparse_1         | True      | True       | False     |      2.64325e-06 |     0.285294  | Optimization terminated successfully |