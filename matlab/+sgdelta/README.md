# +sgdelta package

Reusable MATLAB package for masked perturbation construction and destabilization analysis.

## API overview

### `build_M_from_ss(P, idx)`
Extracts `M = G_rw` from a block-partitioned LTI model `P` with index struct `idx`.

Returns:
- `M` / `Grw`: transfer from disturbance-like channels (`w`) to regulated outputs (`r`),
- optional block struct `G` with `Gyu`, `Gyw`, `Gru`, `Grw`.

### `mask_M(M, rmask, cmask, doShrink)`
Applies binary row/column masks to `M`.

Supports numeric and LTI inputs. Can optionally shrink to active channels.

### `construct_delta_allpass(M, wgrid, opts)`
Builds necessity-side perturbation `Δ(s)` using real all-pass factors.

Also returns:
- peak frequency `wstar`,
- peak singular value `sigmax`,
- diagnostics in `info`.

### `process_masks(M, rmask, wmask, wgrid)`
Convenience wrapper that:
1. masks `M`,
2. constructs `Δ(s)`,
3. returns masked plant and diagnostics.

### `print_info(M, Delta, wstar, sigmax, info)`
Formatted diagnostic printing utility for quick sanity checks.

### `sweep_masks_and_cross_test(Mfull, idx, wgrid, opts)`
Enumerates row/column masks, constructs perturbations, cross-tests all `(M_j, Δ_i)` pairs, and exports CSV reports.

Primary output is a struct with mask definitions, constructed models, destabilization matrix, and export paths.

## Design notes

- Functionality is organized as composable package functions to keep scripts small.
- Most workflows should call these functions through the `sgdelta` namespace.
- Keep this folder MATLAB-only; place future Python logic in a separate top-level Python package directory.

## Maintenance checklist

When adding package functions:
- Include a MATLAB help header (`% FUNCTION_NAME ...`) with inputs/outputs.
- Prefer deterministic outputs and explicit options structs.
- Add the new function to this README.
- Update root `README.md` repository tree if structure changes.
