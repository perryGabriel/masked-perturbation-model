# Realization inspection: toy3_w3_full / baseline_zero / restart 0

- Feasible: `True`
- Objective: `0`

## Numerical matrices

### Q
```text
[[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
```

### Contracted subsystem P=(I-Q)G
```text
[[0.4  0.   0.  ]
 [0.   0.55 0.  ]
 [0.   0.   0.7 ]]
```

## SymPy forms

```text
Q_sympy = MutableDenseMatrix([[Float('0.0', precision=53), Float('0.0', precision=53), Float('0.0', precision=53)], [Float('0.0', precision=53), Float('0.0', precision=53), Float('0.0', precision=53)], [Float('0.0', precision=53), Float('0.0', precision=53), Float('0.0', precision=53)]])
P_sympy = MutableDenseMatrix([[Float('0.40000000000000002', precision=53), Float('0.0', precision=53), Float('0.0', precision=53)], [Float('0.0', precision=53), Float('0.55000000000000004', precision=53), Float('0.0', precision=53)], [Float('0.0', precision=53), Float('0.0', precision=53), Float('0.69999999999999996', precision=53)]])
```
