function out = sweep_masks_and_cross_test(Mfull, idx, wgrid, opts)
% SWEEP_MASKS_AND_CROSS_TEST
% Enumerate binary row/col masks on M = G_rw, build Δ_i for each masked M_i,
% cross-test every Δ_i against every M_j, and export CSVs.
%
% destab(j,i) is TRUE iff either:
%   (A) Den_{j,i}(s) = (I - M_j(s) Δ_i(s))^{-1} has a pole with Re >= stab_eps
% or
%   (B)  min_ω σ_min( I - M_j(jω) Δ_i(jω) ) < ret_tol       (frequency test)
%
% Options (all optional):
%   .min_ones_r, .max_ones_r, .min_ones_c, .max_ones_c   (mask pruning)
%   .stab_eps   : real-part threshold, default -1e-8
%   .ret_tol    : σ_min threshold for frequency test, default 1e-6
%   .sigma_grid : frequency grid (rad/s) for frequency test (default = wgrid)
%   .verbose    : print progress (default true)
%   .csv_prefix : file prefix for exports (default 'mask_sweep')
%
% Requires:
%   - construct_delta_allpass.m  (v2)
%   - close_loop_demo.m / pad_delta.m NOT required here

  if nargin < 3 || isempty(wgrid), wgrid = logspace(-2,2,2001); end
  if nargin < 4, opts = struct; end
  if ~isfield(opts,'min_ones_r'),  opts.min_ones_r = 1; end
  if ~isfield(opts,'max_ones_r'),  opts.max_ones_r = Inf; end
  if ~isfield(opts,'min_ones_c'),  opts.min_ones_c = 1; end
  if ~isfield(opts,'max_ones_c'),  opts.max_ones_c = Inf; end
  if ~isfield(opts,'stab_eps'),    opts.stab_eps   = -1e-8; end
  if ~isfield(opts,'ret_tol'),     opts.ret_tol    = 1e-6; end
  if ~isfield(opts,'sigma_grid'),  opts.sigma_grid = wgrid; end
  if ~isfield(opts,'verbose'),     opts.verbose    = true; end
  if ~isfield(opts,'csv_prefix'),  opts.csv_prefix = 'mask_sweep'; end

  % === Build M = G_rw and sizes ===
  q = size(Mfull,1);         % rows = r-channels
  p = size(Mfull,2);         % cols = w-channels

  % === Enumerate masks ===
  R = enumerate_masks(q, opts.min_ones_r, opts.max_ones_r);
  C = enumerate_masks(p, opts.min_ones_c, opts.max_ones_c);
  Nr = size(R,1); Nc = size(C,1); Nm = Nr*Nc;

  if opts.verbose
    fprintf('Enumerating %d row masks × %d col masks = %d total.\n', Nr, Nc, Nm);
  end

  masks_r = false(Nm, q);
  masks_c = false(Nm, p);
  Mmask   = cell(Nm,1);
  Delta   = cell(Nm,1);
  wstar   = NaN(Nm,1);
  sigmax  = NaN(Nm,1);
  self_destab = false(Nm,1);
  self_min_sigma = NaN(Nm,1);

  % diagonal mask builders (as LTI identity with 0/1 entries)
  mkDiag = @(v) diag(double(v(:)));
  Iq = tf(eye(q));

  % === Construct Δ_i for each masked M_i (no shrinking) ===
  k = 0;
  for ir = 1:Nr
    for ic = 1:Nc
      k = k + 1;
      rmask = R(ir,:).';
      cmask = C(ic,:).';
      masks_r(k,:) = rmask.';  masks_c(k,:) = cmask.';

      Dl = mkDiag(rmask);  Dr = mkDiag(cmask);
      Mi = Dl * Mfull * Dr;             % keep full size; zeroed channels are off

      [Di, w_i, s_i, info_i] = construct_delta_allpass(Mi, wgrid, struct('verbose',false));
      Mmask{k} = Mi;  Delta{k} = Di;  wstar(k) = w_i;  sigmax(k) = s_i;

      % Self check: poles and frequency σmin
      Den_ii = safe_resolvent(Iq, Mi, Di);
      pol_ii = pole(Den_ii);
      is_unstable_poles = any(real(pol_ii) >= opts.stab_eps);
      self_min_sigma(k) = min_sigma_ret(Mi, Di, opts.sigma_grid);
      self_destab(k) = is_unstable_poles || (self_min_sigma(k) < opts.ret_tol);

      if opts.verbose
        fprintf('i=%d/%d  ones(r)=%d ones(c)=%d  w*=%.4g  σ*=%.4g  self: poles=%d  minσ=%.2e\n', ...
          k, Nm, sum(rmask), sum(cmask), w_i, s_i, is_unstable_poles, self_min_sigma(k));
      end
    end
  end

  % === Cross-test every Δ_i vs every M_j ===
  destab = false(Nm, Nm);
  min_sigma = NaN(Nm, Nm);   % min σmin over ω for each (j,i), optional but useful

  if opts.verbose
    fprintf('Enumerating %d cross-combinations.\n', Nm);
  end
  for i = 1:Nm
    Di = Delta{i};
    if opts.verbose
      fprintf('i=%d/%d\n', i, Nm);
    end
    for j = 1:Nm
      Mj = Mmask{j};
      Den_ji = safe_resolvent(Iq, Mj, Di);
      pj = pole(Den_ji);
      bad_poles = any(real(pj) >= opts.stab_eps);
      ms = min_sigma_ret(Mj, Di, opts.sigma_grid);
      min_sigma(j,i) = ms;
      destab(j,i) = bad_poles || (ms < opts.ret_tol);
    end
  end

  % === Export CSVs with BOTH stable & unstable entries ===
  T_masks = table((1:Nm).', masks_r, masks_c, wstar, sigmax, self_destab, self_min_sigma, ...
                  'VariableNames', {'mask_id','rmask','cmask','wstar','sigmax','self_destab','self_min_sigma'});

  [J,I] = ndgrid(1:Nm, 1:Nm);   % J=row=M_j idx, I=col=Δ_i idx
  T_edges = table(I(:), J(:), destab(:), min_sigma(:), ...
                  'VariableNames', {'delta_id','M_id','destabilizes','min_sigma'});

  masks_csv = sprintf('%s_masks.csv', opts.csv_prefix);
  edges_csv = sprintf('%s_destab.csv', opts.csv_prefix);
  writetable(T_masks, masks_csv);
  writetable(T_edges, edges_csv);

  if opts.verbose
    frac = mean(destab(:));
    fprintf('Cross-test matrix: %d of %d pairs destabilized (%.1f%%)\n', ...
      nnz(destab), numel(destab), 100*frac);
    fprintf('Wrote %s and %s\n', masks_csv, edges_csv);
  end

  % === Return struct ===
  out = struct('masks_r', masks_r, 'masks_c', masks_c, ...
               'Delta', {Delta}, 'Mmask', {Mmask}, ...
               'wstar', wstar, 'sigmax', sigmax, ...
               'destab', destab, 'min_sigma', min_sigma, ...
               'self_destab', self_destab, 'self_min_sigma', self_min_sigma, ...
               'masks_csv', masks_csv, 'edges_csv', edges_csv);
end

% ======= helpers =======

function D = safe_resolvent(Iq, M, Delta)
% Den(s) = (I - M(s)Δ(s))^{-1}. Prefer feedback; otherwise inv.
  try
    D = feedback(Iq, M*Delta, +1);    % (I - MΔ)^(-1)
  catch
    D = inv(Iq - M*Delta);            % fallback (can be touchy if improper)
  end
end

function v = min_sigma_ret(M, Delta, wgrid)
% Return min_ω σ_min(I - M(jω) Δ(jω)) over wgrid
  v = inf;
  for k = 1:numel(wgrid)
    jw = 1j*wgrid(k);
    H = eye(size(M,1)) - evalfr(M, jw) * evalfr(Delta, jw);
    s = min(svd(H));
    if s < v, v = s; end
  end
end

function R = enumerate_masks(n, min1, max1)
% All logical masks length n with #ones in [min1, max1]
  if isinf(max1), max1 = n; end
  N = 2^n;
  R = false(0,n);
  for v = 0:N-1
    bits = dec2bin(v, n) == '1';
    k = nnz(bits);
    if k >= min1 && k <= max1
      R(end+1,:) = bits; %#ok<AGROW>
    end
  end
end
