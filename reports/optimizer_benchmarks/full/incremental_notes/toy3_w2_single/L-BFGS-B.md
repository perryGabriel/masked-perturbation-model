# Optimization Notes: toy3_w2_single / L-BFGS-B

- Total runs so far: 4
- Feasible runs so far: 4
- Best feasible objective: 0.901532
- Best restart id: 2
- Best runtime (sec): 0.247

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |     11 | zero             | True      | True       | False     |         0.905068 |      0.233768 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.907043 |      0.329824 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.901532 |      0.246987 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.901532 |      0.227154 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |