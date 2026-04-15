# scripts/

Experiment drivers and utility scripts for the MATLAB implementation.

These scripts are intentionally lightweight wrappers around `+sgdelta` package functions. Keep reusable logic in package code and use scripts for orchestration, plotting, and exports.

## Scripts

### `ACC_Example.m`
Demonstration script that:
- builds a simple transfer-function-based plant/controller configuration,
- constructs multiple masks and perturbations,
- and evaluates closed-loop pole stability for combinations of plants and perturbations.

Use this as a sandbox example, not a production pipeline.

### `four_by_four_ex.m`
Constructs a 4×4 state-space system and demonstrates extraction of `M = G_rw`.

Useful for validating indexing (`idx`) and matrix/channel conventions.

### `get_csv.m`
Runs `sweep_masks_and_cross_test` with configurable options and writes two CSVs:
- `<prefix>_masks.csv`
- `<prefix>_destab.csv`

Use this script when you want a quick export for external analysis.

### `panic_save.m`
Emergency backup utility that snapshots an in-memory `out` struct from a mask sweep run.

Writes:
- MAT archive (highest fidelity),
- tabular CSV summaries,
- optional dense CSV matrices,
- compact JSON metadata.

Run this when a long experiment has completed and you want immediate redundancy before additional processing.

## Conventions for new scripts

- Start each script with a short purpose comment.
- Assume repository root as working directory unless documented otherwise.
- Prefer calls to `sgdelta.<function>` over unqualified function names.
- Keep side effects explicit (file outputs, plots, workspace variables).
- Document generated artifacts near the write calls.
