# Realization inspection: toy2_w1_full / Nelder-Mead / restart 3

- Parameterization: `static_hollow`
- Feasible: `True`
- Objective: `0.702491`

## Numerical matrices

### Q
```text
[[ 0.      -0.21116]
 [-0.19707  0.     ]]
```

### Contracted subsystem P=(I-Q)G
```text
[[0.4     0.14781]
 [0.07883 0.7    ]]
```

## SymPy forms

```text
Q_sympy = MutableDenseMatrix([[Float('0.0', precision=53), Float('-0.21116094042264164', precision=53)], [Float('-0.19707125442497822', precision=53), Float('0.0', precision=53)]])
P_sympy = MutableDenseMatrix([[Float('0.40000000000000002', precision=53), Float('0.14781265829584914', precision=53)], [Float('0.078828501769991288', precision=53), Float('0.69999999999999996', precision=53)]])
```
