# Optimization Notes: toy2_w1_full / basinhopping

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 0.702491
- Best restart id: 0
- Best runtime (sec): 0.5147

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                                                |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.702491 |      0.514701 | ['requested number of basinhopping iterations completed successfully'] |
|            1 |     24 | gaussian_0       | False     | False      | True      |       inf        |      5        | timeout                                                                |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.702491 |      0.391656 | ['requested number of basinhopping iterations completed successfully'] |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.702491 |      0.414916 | ['requested number of basinhopping iterations completed successfully'] |