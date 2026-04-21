# Optimization Notes: toy3_w2_single / L-BFGS-B

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 0.901532
- Best restart id: 2
- Best runtime (sec): 0.2049

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |       inf        |      5        | timeout                                              |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.907043 |      0.328218 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.901532 |      0.204925 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.901532 |      0.191637 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |