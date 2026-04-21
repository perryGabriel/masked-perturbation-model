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

## v1 assumptions for `Q`/`P` parameterization

For the current implementation generation (v1), we use a constrained `Q` parameterization:

- For `i != j`, each `Q_{ij}(z)` is modeled as a discrete-time proper rational transfer function.
- Diagonal entries are fixed to zero exactly: `Q_{ii}(z) = 0`.
- A shared denominator order `d_Q` is used across all nonzero off-diagonal entries (recommended v1 default).
- v1 chooses **strictly proper** off-diagonal entries:
  - `deg(num(Q_{ij})) <= d_Q - 1` for all modeled links.
- `P` is never parameterized independently; it is always derived/computed from `Q` and `G`.

Tradeoff (why shared denominator in v1):
- reduces the decision-variable count;
- improves practical identifiability by limiting equivalent parameterizations;
- tends to improve conditioning in numerical optimization and estimation routines.
