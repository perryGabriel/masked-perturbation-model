# Optimization Notes: toy3_w2_single / L-BFGS-B

- Total runs so far: 8
- Feasible runs so far: 8
- Best feasible objective: 0.897716
- Best restart id: 5
- Best runtime (sec): 0.3344

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                              |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:-----------------------------------------------------|
|            0 |    202 | zero             | True      | True       | False     |         0.905068 |      0.390016 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            1 |    215 | gaussian_0       | True      | True       | False     |         0.904445 |      0.301795 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            2 |    228 | sparse_0         | True      | True       | False     |         0.904445 |      1.11468  | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            3 |    241 | gaussian_1       | True      | True       | False     |         0.901826 |      0.412297 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            4 |    254 | sparse_1         | True      | True       | False     |         0.901826 |      0.238246 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            5 |    267 | gaussian_2       | False     | True       | False     |         0.897716 |      0.334423 | ABNORMAL:                                            |
|            6 |    280 | sparse_2         | True      | True       | False     |         0.90833  |      0.491892 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
|            7 |    293 | gaussian_3       | True      | True       | False     |         0.909242 |      0.831369 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |