# Optimization Notes: toy3_w2_single / L-BFGS-B

- Total runs so far: 5
- Feasible runs so far: 5
- Best feasible objective: 0.897716
- Best restart id: 4
- Best runtime (sec): 0.4247

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |      7 | zero             | True      | True       | False     |         0.905068 |      0.484021 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            1 |     20 | gaussian_0       | True      | True       | False     |         0.907043 |      0.37244  | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |     33 | sparse_0         | True      | True       | False     |         0.907043 |      0.567206 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            3 |     46 | gaussian_1       | True      | True       | False     |         0.897716 |      0.453377 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL     |
|            4 |     59 | sparse_1         | True      | True       | False     |         0.897716 |      0.424683 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |