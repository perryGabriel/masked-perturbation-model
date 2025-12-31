% Assume you already have P, idx built, plus construct_delta_allpass.m on path.
wgrid = logspace(-2,2,2001);
opts = struct('min_ones_r', 1, 'max_ones_r', 2, ...
              'min_ones_c', 1, 'max_ones_c', 2, ...
              'stab_eps', -1e-8, 'verbose', true, ...
              'csv_prefix', 'mask_sweep_run1');

out = sweep_masks_and_cross_test(M, idx, wgrid, opts);

% CSVs:
%   - 'mask_sweep_run1_masks.csv'   (one row per mask, with rmask/cmask cells)
%   - 'mask_sweep_run1_destab.csv'  (rows: delta_id, M_id, destabilizes=true)
