# Optimization Notes: toy2_w1_full / Nelder-Mead

- Total runs so far: 4
- Feasible runs so far: 2
- Best feasible objective: 0.702491
- Best restart id: 3
- Best runtime (sec): 0.01659

## Failure modes observed

timeout: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                               |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |       inf        |     5         | timeout                               |
|            1 |     24 | gaussian_0       | False     | False      | True      |       inf        |     5         | timeout                               |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.702491 |     0.0137255 | Optimization terminated successfully. |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.702491 |     0.0165912 | Optimization terminated successfully. |