# Chapter-2 feedback vulnerability (FMP) extension

This repository now includes a chapter-2 numerical exploration toolkit for small systems.

## Scope and honesty

- These utilities provide **numerical sanity checks** on toy examples.
- They do **not** prove theorems in full generality.
- The design methods are explicitly finite-dimensional surrogates/relaxations.

## Measured vulnerability

For a chosen system and realization matrix `Q`, we evaluate:

- Full attacker: `V_full(Q) = ||M_Q||` (spectral norm for static demos)
- Single-link attacker: `V_single(Q) = max_ij |(M_Q)_{ij}|`

with

`M_Q = Cbar * Gamma * (I - Q)^(-1) * Pbar`,
`Gamma = (I - G alpha)^(-1)`.

## Bounds implemented

- `bound_upper_full` (Q=0 helper used as an experiment-side overlay).
- `bound_lower_single_link` (under chapter assumptions; static simplification).
- `bound_lower_full` (static simplification using `sigma_min(Pbar)/(1+||G alpha||)`).

The code avoids silent cross-threat comparison by labeling results and plotting per threat model.

## Design methods

- `projection_design`: iterative regularized least-squares surrogate inspired by projection derivation.
- `lp_relaxation_design`: LP-inspired static linearization/relaxation using `scipy.optimize.linprog`.

Both methods report true measured vulnerability after design and should be interpreted as approximations.

## v1 modeling assumptions for `Q`

To keep chapter-2 experiments well-posed and numerically stable, v1 uses the following structural assumptions:

- **Discrete-time rational parameterization:** for each off-diagonal link `i != j`, `Q_{ij}(z)` is a **proper** discrete-time rational transfer function.
- **No self-loops in `Q`:** `Q_{ii}(z) = 0` exactly for all `i`.
- **Shared denominator order (recommended and used in v1):** all nonzero off-diagonal `Q_{ij}(z)` share a common denominator polynomial order `d_Q`.
- **Numerator degree choice for v1:** we enforce **strictly proper** entries, i.e., `deg(num(Q_{ij})) <= d_Q - 1`.
- **No independent `P` parameterization:** `P` is always computed from `Q` and `G` (not optimized as an independent decision variable).

### Why the shared denominator is recommended in v1

Using one denominator basis across off-diagonal entries reduces the number of free coefficients and typically improves numerical behavior:

- fewer decision variables in optimization;
- less over-parameterization, improving identifiability;
- better conditioning for regression/projection-style updates.

This is a pragmatic v1 choice to improve robustness of small/medium-scale numerical experiments.

## How to run

- Notebook bounds demo: `notebooks/chapter2_bounds_demo.ipynb`
- Notebook design demo: `notebooks/chapter2_design_demo.ipynb`
- Tests: `pytest`

## DSF migration note

- DSF usage is opt-in via explicit names and constructors (for example `benchmark_problem_registry_dsf(...)` / `_toy_problem_dsf(...)`).
- The default benchmark registry continues to use static hollow parameterization for backward-compatible reproductions.
- “No behavior change unless `ParameterizationSpec.kind='dsf_poly'`”.
