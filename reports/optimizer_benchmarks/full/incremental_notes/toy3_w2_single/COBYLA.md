# Optimization Notes: toy3_w2_single / COBYLA

- Total runs so far: 4
- Feasible runs so far: 4
- Best feasible objective: 0.92207
- Best restart id: 3
- Best runtime (sec): 0.05043

## Failure modes observed

No infeasible runs for this algorithm/problem slice.

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                                   |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:----------------------------------------------------------|
|            0 |     11 | zero             | False     | True       | False     |         0.926966 |     0.0602242 | Maximum number of function evaluations has been exceeded. |
|            1 |     24 | gaussian_0       | True      | True       | False     |         0.954904 |     0.0654074 | Optimization terminated successfully.                     |
|            2 |     37 | sparse_0         | True      | True       | False     |         0.941055 |     0.0532826 | Optimization terminated successfully.                     |
|            3 |     50 | gaussian_1       | True      | True       | False     |         0.92207  |     0.0504303 | Optimization terminated successfully.                     |