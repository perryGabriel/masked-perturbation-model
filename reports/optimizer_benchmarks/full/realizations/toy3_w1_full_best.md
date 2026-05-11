# Realization inspection: toy3_w1_full / differential_evolution / restart 1

- Parameterization: `static_hollow`
- Feasible: `True`
- Objective: `0.706888`

## Numerical matrices

### Q
```text
[[ 0.       0.2486  -0.03809]
 [ 0.08979  0.      -0.04886]
 [ 0.0655  -0.15426  0.     ]]
```

### Contracted subsystem P=(I-Q)G
```text
[[ 0.4     -0.13673  0.02667]
 [-0.03591  0.55     0.0342 ]
 [-0.0262   0.08484  0.7    ]]
```

## SymPy forms

```text
Q_sympy = MutableDenseMatrix([[Float('0.0', precision=53), Float('0.24860120890407073', precision=53), Float('-0.038093140649819807', precision=53)], [Float('0.089785949982013835', precision=53), Float('0.0', precision=53), Float('-0.048858217635723566', precision=53)], [Float('0.065496732210704844', precision=53), Float('-0.15425754405171721', precision=53), Float('0.0', precision=53)]])
P_sympy = MutableDenseMatrix([[Float('0.40000000000000002', precision=53), Float('-0.1367306648972389', precision=53), Float('0.026665198454873862', precision=53)], [Float('-0.035914379992805535', precision=53), Float('0.55000000000000004', precision=53), Float('0.034200752345006497', precision=53)], [Float('-0.026198692884281938', precision=53), Float('0.084841649228444477', precision=53), Float('0.69999999999999996', precision=53)]])
```
