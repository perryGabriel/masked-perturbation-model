# Realization inspection: toy3_w2_single / L-BFGS-B / restart 4

- Feasible: `True`
- Objective: `0.897716`

## Numerical matrices

### Q
```text
[[ 0.       0.25     0.25   ]
 [-0.17644  0.       0.25   ]
 [-0.25    -0.25     0.     ]]
```

### Contracted subsystem P=(I-Q)G
```text
[[ 0.4     -0.1375  -0.175  ]
 [ 0.07058  0.55    -0.175  ]
 [ 0.1      0.1375   0.7    ]]
```

## SymPy forms

```text
Q_sympy = MutableDenseMatrix([[Float('0.0', precision=53), Float('0.25', precision=53), Float('0.25', precision=53)], [Float('-0.17644141671988844', precision=53), Float('0.0', precision=53), Float('0.25', precision=53)], [Float('-0.25', precision=53), Float('-0.25', precision=53), Float('0.0', precision=53)]])
P_sympy = MutableDenseMatrix([[Float('0.40000000000000002', precision=53), Float('-0.13750000000000001', precision=53), Float('-0.17499999999999999', precision=53)], [Float('0.07057656668795538', precision=53), Float('0.55000000000000004', precision=53), Float('-0.17499999999999999', precision=53)], [Float('0.10000000000000001', precision=53), Float('0.13750000000000001', precision=53), Float('0.69999999999999996', precision=53)]])
```
