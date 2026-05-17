# Optimization Notes: toy2_w1_full / L-BFGS-B

- Total runs so far: 2
- Feasible runs so far: 2
- Best feasible objective: 0.702491
- Best restart id: 0
- Best runtime (sec): 0.005631

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                          |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-------------------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.702491 |     0.0056308 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.702491 |     0.0239523 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |