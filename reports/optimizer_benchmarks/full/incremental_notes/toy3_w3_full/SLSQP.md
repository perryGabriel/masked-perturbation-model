# Optimization Notes: toy3_w3_full / SLSQP

- Total runs so far: 4
- Feasible runs so far: 2
- Best feasible objective: 1.90081e-06
- Best restart id: 1
- Best runtime (sec): 0.1713

## Failure modes observed

timeout: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |    inf           |      5        | timeout                              |
|            1 |     24 | gaussian_0       | True      | True       | False     |      1.90081e-06 |      0.171324 | Optimization terminated successfully |
|            2 |     37 | sparse_0         | False     | False      | True      |    inf           |      5        | timeout                              |
|            3 |     50 | gaussian_1       | True      | True       | False     |      2.33277e-06 |      0.18316  | Optimization terminated successfully |