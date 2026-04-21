# Optimization Notes: toy3_w2_single / COBYLA

- Total runs so far: 8
- Feasible runs so far: 6
- Best feasible objective: 0.920429
- Best restart id: 5
- Best runtime (sec): 0.03665

## Failure modes observed

Optimization terminated successfully.: 2

## Recent runs

|   restart_id |   seed | initialization   | success   | feasible   | timeout   |   true_objective |   runtime_sec | message                               |
|-------------:|-------:|:-----------------|:----------|:-----------|:----------|-----------------:|--------------:|:--------------------------------------|
|            0 |    202 | zero             | True      | True       | False     |         0.926759 |     0.0293401 | Optimization terminated successfully. |
|            1 |    215 | gaussian_0       | True      | True       | False     |         0.941264 |     0.0178762 | Optimization terminated successfully. |
|            2 |    228 | sparse_0         | True      | False      | False     |         0.911647 |     0.0414631 | Optimization terminated successfully. |
|            3 |    241 | gaussian_1       | True      | False      | False     |         0.960368 |     0.0293282 | Optimization terminated successfully. |
|            4 |    254 | sparse_1         | True      | True       | False     |         0.965686 |     0.0243529 | Optimization terminated successfully. |
|            5 |    267 | gaussian_2       | True      | True       | False     |         0.920429 |     0.036648  | Optimization terminated successfully. |
|            6 |    280 | sparse_2         | True      | True       | False     |         0.942012 |     0.0269384 | Optimization terminated successfully. |
|            7 |    293 | gaussian_3       | True      | True       | False     |         0.956961 |     0.0467828 | Optimization terminated successfully. |