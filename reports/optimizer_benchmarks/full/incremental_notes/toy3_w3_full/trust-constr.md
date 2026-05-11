# Optimization Notes: toy3_w3_full / trust-constr

- Total runs so far: 4
- Feasible runs so far: 1
- Best feasible objective: 7.08452e-10
- Best restart id: 0
- Best runtime (sec): 1.23

## Failure modes observed

timeout: 3

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                                 |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------------------------|
|            0 |     11 | zero             | False     | True       | False     |      7.08452e-10 |        1.2297 | The maximum number of function evaluations is exceeded. |
|            1 |     24 | gaussian_0       | False     | False      | True      |    inf           |        5      | timeout                                                 |
|            2 |     37 | sparse_0         | False     | False      | True      |    inf           |        5      | timeout                                                 |
|            3 |     50 | gaussian_1       | False     | False      | True      |    inf           |        5      | timeout                                                 |