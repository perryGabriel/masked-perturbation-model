# Optimization Notes: toy2_w1_full / SLSQP

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 0.702491
- Best restart id: 1
- Best runtime (sec): 0.009107

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.702491 |     0.004565  | Optimization terminated successfully |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0091068 | Optimization terminated successfully |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.702491 |     0.0078837 | Optimization terminated successfully |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |     5         | timeout                              |