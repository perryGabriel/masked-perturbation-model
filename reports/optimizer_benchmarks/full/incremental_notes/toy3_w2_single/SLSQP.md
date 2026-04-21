# Optimization Notes: toy3_w2_single / SLSQP

- Total runs so far: 4
- Feasible runs so far: 1
- Best feasible objective: 0.905068
- Best restart id: 0
- Best runtime (sec): 0.1173

## Failure modes observed

timeout: 3

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.905068 |      0.117337 | Optimization terminated successfully |
|            1 |     24 | gaussian_0       | False     | False      | True      |       inf        |      5        | timeout                              |
|            2 |     37 | sparse_0         | False     | False      | True      |       inf        |      5        | timeout                              |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |      5        | timeout                              |