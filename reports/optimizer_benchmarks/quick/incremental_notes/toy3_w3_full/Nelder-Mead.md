# Optimization Notes: toy3_w3_full / Nelder-Mead

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0
- Best restart id: 0
- Best runtime (sec): 0.01065

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                         |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:------------------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0        |     0.0106486 | Optimization terminated successfully.           |
|            1 |     20 | gaussian_0       | False     | True       | False     |         0.110197 |     0.0700014 | Maximum number of iterations has been exceeded. |
|            2 |     33 | sparse_0         | False     | True       | False     |         0.001493 |     0.0894129 | Maximum number of iterations has been exceeded. |
|            3 |     46 | gaussian_1       | False     | True       | False     |         0.113127 |     0.0829273 | Maximum number of iterations has been exceeded. |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.012349 |     0.0819017 | Optimization terminated successfully.           |