# Optimization Notes: toy2_w1_full / COBYLA

- Total runs so far: 4
- Feasible runs so far: 1
- Best feasible objective: 0.702491
- Best restart id: 1
- Best runtime (sec): 0.03261

## Failure modes observed

timeout: 3

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                               |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |       inf        |      5        | timeout                               |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.702491 |      0.032606 | Optimization terminated successfully. |
|            2 |     37 | sparse_0         | False     | False      | True      |       inf        |      5        | timeout                               |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |      5        | timeout                               |