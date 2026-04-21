# Realization inspection: toy3_w1_full / Nelder-Mead / restart 7

- Feasible: `True`
- Objective: `0.706888`

## Numerical matrices

### Q
```text
[[ 0.      -0.23562  0.0725 ]
 [ 0.23965  0.      -0.22694]
 [ 0.20439 -0.17259  0.     ]]
```

### Contracted subsystem P=(I-Q)G
```text
[[ 0.4      0.12959 -0.05075]
 [-0.09586  0.55     0.15886]
 [-0.08176  0.09492  0.7    ]]
```

## SymPy forms

```text
Q_sympy = MutableDenseMatrix([[Float('0.0', precision=53), Float('-0.23562124519477826', precision=53), Float('0.072504151866787228', precision=53)], [Float('0.23965134867279891', precision=53), Float('0.0', precision=53), Float('-0.22693723754106948', precision=53)], [Float('0.20439109098311706', precision=53), Float('-0.17258890320331263', precision=53), Float('0.0', precision=53)]])
P_sympy = MutableDenseMatrix([[Float('0.40000000000000002', precision=53), Float('0.12959168485712805', precision=53), Float('-0.050752906306751056', precision=53)], [Float('-0.095860539469119571', precision=53), Float('0.55000000000000004', precision=53), Float('0.15885606627874863', precision=53)], [Float('-0.081756436393246826', precision=53), Float('0.09492389676182196', precision=53), Float('0.69999999999999996', precision=53)]])
```
