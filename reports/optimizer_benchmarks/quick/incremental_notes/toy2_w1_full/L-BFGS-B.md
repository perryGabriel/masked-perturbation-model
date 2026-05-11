# Optimization Notes: toy2_w1_full / L-BFGS-B

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.702491
- Best restart id: 3
- Best runtime (sec): 0.01305

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                          |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.702491 |     0.0093243 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0058344 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            2 |     33 | sparse_0         | True      | True       | False     |         0.702491 |     0.0093058 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            3 |     46 | gaussian_1       | True      | True       | False     |         0.702491 |     0.0130476 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.702491 |     0.0090867 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |