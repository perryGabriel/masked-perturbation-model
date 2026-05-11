# Optimization Notes: toy3_w2_single / dual_annealing

- Total runs so far: 8
- Feasible runs so far: 7
- Best feasible objective: 0.897716
- Best restart id: 5
- Best runtime (sec): 2.677

## Failure modes observed

timeout: 1

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                                 |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:----------------------------------------|
|            0 |    202 | zero             | True      | True       | False     |         0.905127 |       3.73483 | ['Maximum number of iteration reached'] |
|            1 |    215 | gaussian_0       | True      | True       | False     |         0.897716 |       4.89456 | ['Maximum number of iteration reached'] |
|            2 |    228 | sparse_0         | True      | True       | False     |         0.90648  |       3.6259  | ['Maximum number of iteration reached'] |
|            3 |    241 | gaussian_1       | True      | True       | False     |         0.897716 |       4.99373 | ['Maximum number of iteration reached'] |
|            4 |    254 | sparse_1         | True      | True       | False     |         0.902131 |       4.24866 | ['Maximum number of iteration reached'] |
|            5 |    267 | gaussian_2       | True      | True       | False     |         0.897716 |       2.67749 | ['Maximum number of iteration reached'] |
|            6 |    280 | sparse_2         | True      | True       | False     |         0.90833  |       5.30515 | ['Maximum number of iteration reached'] |
|            7 |    293 | gaussian_3       | False     | False      | True      |       inf        |      45       | timeout                                 |