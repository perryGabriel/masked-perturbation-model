# Optimization Notes: toy3_w3_full / L-BFGS-B

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0
- Best restart id: 0
- Best runtime (sec): 0.1768

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |      7 | zero             | False     | True       | False     |      0           |      0.176782 | ABNORMAL:                                            |
|            1 |     20 | gaussian_0       | True      | True       | False     |      2.32783e-08 |      1.04122  | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |     33 | sparse_0         | False     | True       | False     |      6.47545e-09 |      0.350774 | ABNORMAL:                                            |
|            3 |     46 | gaussian_1       | True      | True       | False     |      1.12538e-08 |      0.507277 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            4 |     59 | sparse_1         | True      | True       | False     |      8.20817e-09 |      0.647887 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |