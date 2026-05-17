# Optimization Notes: toy3_w1_full / trust-constr

- Total runs so far: 10
- Feasible runs so far: 9
- Best feasible objective: 0.706888
- Best restart id: 7
- Best runtime (sec): 0.03035

## Failure modes observed

BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                                                                                           |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:------------------------------------------------------------------------------------------------------------------|
|            0 |    101 | zero             | False     | False      | False     |       inf        |     0         | BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending. |
|            1 |    114 | gaussian_0       | True      | True       | False     |         0.706888 |     0.0495432 | `gtol` termination condition is satisfied.                                                                        |
|            2 |    127 | sparse_0         | True      | True       | False     |         0.706888 |     0.0273732 | `gtol` termination condition is satisfied.                                                                        |
|            3 |    140 | gaussian_1       | True      | True       | False     |         0.706888 |     0.0227861 | `gtol` termination condition is satisfied.                                                                        |
|            4 |    153 | sparse_1         | True      | True       | False     |         0.706888 |     0.0145248 | `gtol` termination condition is satisfied.                                                                        |
|            5 |    166 | gaussian_2       | True      | True       | False     |         0.706888 |     0.022501  | `gtol` termination condition is satisfied.                                                                        |
|            6 |    179 | sparse_2         | True      | True       | False     |         0.706888 |     0.0296496 | `gtol` termination condition is satisfied.                                                                        |
|            7 |    192 | gaussian_3       | True      | True       | False     |         0.706888 |     0.0303481 | `gtol` termination condition is satisfied.                                                                        |
|            8 |    205 | sparse_3         | True      | True       | False     |         0.706888 |     0.0229404 | `gtol` termination condition is satisfied.                                                                        |
|            9 |    218 | gaussian_4       | True      | True       | False     |         0.706888 |     0.0210196 | `gtol` termination condition is satisfied.                                                                        |