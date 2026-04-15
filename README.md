# masked-perturbation-model

MATLAB reference implementation for masked perturbation analysis and destabilization cross-testing.

This repository currently contains:
- A namespaced MATLAB package (`+sgdelta`) with reusable core functions.
- Example and utility scripts (`scripts/`) for experiments and data export.
- A lightweight setup helper (`setup.m`) to configure the MATLAB path.

The project is being prepared for a future Python package extension. This README documents the current MATLAB layout so the codebase remains organized as new language-specific modules are added.

## Repository layout

```text
masked-perturbation-model/
├── +sgdelta/                    # Reusable MATLAB package functions
│   ├── build_M_from_ss.m
│   ├── construct_delta_allpass.m
│   ├── mask_M.m
│   ├── print_info.m
│   ├── process_masks.m
│   └── sweep_masks_and_cross_test.m
├── scripts/                     # Reproducible examples and one-off utilities
│   ├── ACC_Example.m
│   ├── four_by_four_ex.m
│   ├── get_csv.m
│   └── panic_save.m
├── setup.m                      # Adds project folders to MATLAB path
└── README.md
```

## MATLAB prerequisites

- MATLAB R2020b+ recommended.
- Control System Toolbox required (`tf`, `ss`, `feedback`, `evalfr`, etc.).

## Quick start

1. Open MATLAB at the repository root.
2. Run setup once per session:

```matlab
setup
```

3. Call package functions using the `sgdelta` namespace:

```matlab
wgrid = logspace(-2, 2, 2001);
[M_d, Delta, wstar, sigmax, info] = sgdelta.process_masks(M, rmask, wmask, wgrid);
```

> **Note:** Legacy scripts may call package functions without the namespace. Prefer `sgdelta.<function>` for new work to avoid name collisions as the repository grows.

## Core workflow (MATLAB)

Typical workflow for one plant model:

1. Build or load an LTI model `P` and index partitions with `idx`.
2. Extract `M = G_rw` using `sgdelta.build_M_from_ss`.
3. Apply row/column masks with `sgdelta.mask_M` (or directly via `sgdelta.process_masks`).
4. Construct necessity-side perturbation `Δ(s)` via `sgdelta.construct_delta_allpass`.
5. Run exhaustive mask/perturbation cross-tests with `sgdelta.sweep_masks_and_cross_test`.
6. Export CSV summaries for downstream analysis.

## Script usage

See `scripts/README.md` for script-by-script descriptions, expected inputs, and outputs.

## Function reference

See `+sgdelta/README.md` for the package API and function-level guidance.

## Project conventions (for upcoming Python extension)

To keep MATLAB and Python code cleanly separated:

- Keep reusable MATLAB functions in `+sgdelta/`.
- Keep MATLAB scripts in `scripts/` and treat them as experiment drivers.
- Avoid adding nontrivial logic directly to top-level scripts; prefer package functions.
- Use explicit namespaces (`sgdelta.<name>`) in new MATLAB code.
- Add new Python code under a dedicated Python directory (e.g., `python/` or `src/`) rather than mixing it into MATLAB folders.

## Documentation quality checklist

Before adding new modules:

- [ ] New files include a short header comment describing purpose and inputs/outputs.
- [ ] Any new script is listed in `scripts/README.md`.
- [ ] Any new package function is listed in `+sgdelta/README.md`.
- [ ] Root README layout stays up to date.
- [ ] Example commands remain runnable from repository root.

## License

No license file is currently included. Add a project license before wider distribution.
