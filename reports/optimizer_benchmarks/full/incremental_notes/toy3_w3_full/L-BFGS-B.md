# Optimization Notes: toy3_w3_full / L-BFGS-B

- Total runs so far: 4
- Feasible runs so far: 3
- Best feasible objective: 9.80942e-09
- Best restart id: 2
- Best runtime (sec): 0.5115

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |     11 | zero             | False     | False      | True      |    inf           |      5        | timeout                                              |
|            1 |     24 | gaussian_0       | True      | True       | False     |      2.71088e-08 |      0.460326 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |     37 | sparse_0         | True      | True       | False     |      9.80942e-09 |      0.511517 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            3 |     50 | gaussian_1       | True      | True       | False     |      2.56299e-08 |      0.446343 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |