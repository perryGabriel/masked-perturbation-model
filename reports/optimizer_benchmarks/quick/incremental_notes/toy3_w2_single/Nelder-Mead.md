# Optimization Notes: toy3_w2_single / Nelder-Mead

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.947919
- Best restart id: 3
- Best runtime (sec): 0.07888

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                         |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:------------------------------------------------|
|            0 |      7 | zero             | False     | True       | False     |         1.00184  |     0.0570369 | Maximum number of iterations has been exceeded. |
|            1 |     20 | gaussian_0       | False     | True       | False     |         0.962878 |     0.0698481 | Maximum number of iterations has been exceeded. |
|            2 |     33 | sparse_0         | True      | True       | False     |         1.00119  |     0.0656513 | Optimization terminated successfully.           |
|            3 |     46 | gaussian_1       | False     | True       | False     |         0.947919 |     0.0788753 | Maximum number of iterations has been exceeded. |
|            4 |     59 | sparse_1         | True      | True       | False     |         1.00014  |     0.0669754 | Optimization terminated successfully.           |