# Optimization Notes: toy2_w1_full / L-BFGS-B

- Total runs so far: 4
- Feasible runs so far: 4
- Best feasible objective: 0.702491
- Best restart id: 1
- Best runtime (sec): 0.002207

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                          |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.702491 |     0.0019364 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0022065 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.702491 |     0.0016546 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.702491 |     0.0016738 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |