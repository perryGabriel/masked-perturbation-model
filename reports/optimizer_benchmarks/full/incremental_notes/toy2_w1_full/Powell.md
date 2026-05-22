# Optimization Notes: toy2_w1_full / Powell

- Total runs so far: 4
- Feasible runs so far: 2
- Best feasible objective: 0.702491
- Best restart id: 0
- Best runtime (sec): 0.03109

## Failure modes observed

timeout: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                               |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.702491 |     0.0310932 | Optimization terminated successfully. |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.702491 |     0.015375  | Optimization terminated successfully. |
|            2 |     37 | sparse_0         | False     | False      | True      |       inf        |     5         | timeout                               |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |     5         | timeout                               |