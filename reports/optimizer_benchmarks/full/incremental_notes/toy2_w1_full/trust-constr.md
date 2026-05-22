# Optimization Notes: toy2_w1_full / trust-constr

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 0.702491
- Best restart id: 1
- Best runtime (sec): 0.008123

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                    |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.702491 |     0.0362993 | `gtol` termination condition is satisfied. |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0081231 | `gtol` termination condition is satisfied. |
|            2 |     37 | sparse_0         | False     | False      | True      |       inf        |     5         | timeout                                    |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.702491 |     0.0074982 | `gtol` termination condition is satisfied. |