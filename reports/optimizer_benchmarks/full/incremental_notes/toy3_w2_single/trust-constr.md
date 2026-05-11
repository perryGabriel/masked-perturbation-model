# Optimization Notes: toy3_w2_single / trust-constr

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 0.903246
- Best restart id: 0
- Best runtime (sec): 1.507

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                                 |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------------------------|
|            0 |     11 | zero             | False     | True       | False     |         0.903246 |       1.50662 | The maximum number of function evaluations is exceeded. |
|            1 |     24 | gaussian_0       | False     | True       | False     |         0.930072 |       1.23486 | The maximum number of function evaluations is exceeded. |
|            2 |     37 | sparse_0         | False     | True       | False     |         0.912962 |       1.0051  | The maximum number of function evaluations is exceeded. |
|            3 |     50 | gaussian_1       | False     | False      | True      |       inf        |       5       | timeout                                                 |