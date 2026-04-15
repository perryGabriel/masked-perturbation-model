# Implementation notes

`mpmgame` is intentionally transfer-function-centric.

- Public APIs accept transfer matrices `M`, attacks `Δ`, and binary masks `∇`.
- Stability checks use poles of `(I - M(Δ∘∇))^{-1}`.
- The alternative masked-model perspective uses:
  - `M_∇ = diag(η_r) M diag(η_w)`
  - with rank-1 mask `∇ = η_w η_r^T`
- The package avoids μ-analysis and focuses on finite action sets, success sets, strategy dominance, and zero-sum equilibria.

Numerical caveat:
- Near-marginal poles trigger warnings for reproducibility and debugging.
