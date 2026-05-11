# Optimization Notes: toy3_w1_full / trust-constr

- Total runs so far: 4
- Feasible runs so far: 2
- Best feasible objective: 0.706888
- Best restart id: 0
- Best runtime (sec): 0.03815

## Failure modes observed

timeout: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                    |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.706888 |     0.0381526 | `gtol` termination condition is satisfied. |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.706888 |     0.0087344 | `gtol` termination condition is satisfied. |
|            2 |     37 | sparse_0         | False     | False      | True      |       inf        |     5         | timeout                                    |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |     5         | timeout                                    |