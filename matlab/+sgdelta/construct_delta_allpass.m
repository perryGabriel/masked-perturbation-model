function [Delta, wstar, sigmax, info] = construct_delta_allpass(M, wgrid, opts)
% CONSTRUCT_DELTA_ALLPASS (v2)  Necessity-side Δ(s) via real all-pass factors.
% - robust to name shadowing (uses evalfr/FRD indexing)
% - no repmat; vectorized phase handling
% - correct phase map: alpha = -w*cot(phi/2)

  if nargin < 2 || isempty(wgrid), wgrid = logspace(-3,3,2001); end
  if nargin < 3, opts = struct; end
  if ~isfield(opts,'verbose'),   opts.verbose   = false; end
  if ~isfield(opts,'refine'),    opts.refine    = true;  end
  if ~isfield(opts,'refine_bw'), opts.refine_bw = 1.15;  end

  fr = @(sys,w) local_fr(sys,w);           % numeric M(jw)
  ap = @(a_,s) local_allpass(a_,s);        % (α - s)/(α + s)
  s  = tf('s');

  % ---- find peak σ ----
  sigmax = -inf; wstar = NaN; Mstar = []; u = []; v = [];
  for w = wgrid(:).'
    A = fr(M,w);
    [U,S,V] = svd(A,'econ'); s1 = S(1,1);
    if s1 > sigmax, sigmax=s1; wstar=w; Mstar=A; u=U(:,1); v=V(:,1); end
  end
  if opts.refine
    wb = max(wgrid(1), wstar/opts.refine_bw);
    wt = min(wgrid(end), wstar*opts.refine_bw);
    for w = linspace(wb,wt,201)
      A = fr(M,w);
      [U,S,V] = svd(A,'econ'); s1 = S(1,1);
      if s1 > sigmax, sigmax=s1; wstar=w; Mstar=A; u=U(:,1); v=V(:,1); end
    end
  end

  % ---- build Δ(s) = (1/σ) * vecV * vecU ----
  [p,q] = size_local(M);         % p = rows (=|r|), q = cols (=|w|)
  theta_u = angle(u);            % phases of u
  phi_v   = angle(v);

  [phi_u_mod, flip_u] = phase_to_neg_half_vec(-theta_u);
  [phi_v_mod, flip_v] = phase_to_neg_half_vec( phi_v);

  alpha = zeros(numel(u),1);
  beta  = zeros(numel(v),1);

  % Preallocate without repmat
  vecV = tf(zeros(q,1));         % q-by-1 tf array
  for j = 1:q
    beta(j) = phase_to_alpha(phi_v_mod(j), wstar);
    Aj = ap(beta(j), s); if flip_v(j), Aj = -Aj; end
    vecV(j,1) = abs(v(j)) * Aj;
  end

  vecU = tf(zeros(1,p));         % 1-by-p tf array
  for i = 1:p
    alpha(i) = phase_to_alpha(phi_u_mod(i), wstar);
    Ai = ap(alpha(i), s); if flip_u(i), Ai = -Ai; end
    vecU(1,i) = abs(u(i)) * Ai;
  end

  Delta  = (1/sigmax) * (vecV * vecU);      % q-by-p

  % ---- diagnostics ----
  Hstar  = fr(Delta, wstar);
  target = (1/sigmax) * (v * (u'));         % u' = conj-transpose
  rd     = eye(size(Mstar,1)) - Mstar * Hstar;

  info = struct();
  info.u=u; info.v=v; info.Mstar=Mstar;
  info.Delta0 = (v*(u'))/sigmax;
  info.alpha=alpha; info.beta=beta;
  info.phase_u=theta_u; info.phase_v=phi_v;
  info.Hstar=Hstar; info.target=target;
  info.norm_err = norm(Hstar - target, 2);
  info.rd_sigma_min = min(svd(rd));

  if opts.verbose
    fprintf('||M||_inf ~= %.6g at w* = %.6g rad/s\n', sigmax, wstar);
    fprintf('‖Δ(jw*) - (1/σ) v u*‖₂ = %.3e\n', info.norm_err);
    fprintf('σ_min(I - M(jw*)Δ(jw*)) = %.3e\n', info.rd_sigma_min);
  end
end

% ---------- helpers ----------
function H = local_fr(sys, w)
  if isnumeric(sys)
    H = sys;
  elseif isa(sys,'frd')
    [~,k] = min(abs(sys.Frequency - w));
    H = sys.ResponseData(:,:,k);
  else
    H = evalfr(sys, 1j*w);
  end
end

function A = local_allpass(alpha, s)
  TOL = 1e-12;
  if isinf(alpha), A = 1;
  elseif alpha < TOL, A = -1;
  else, A = tf([-1, alpha], [1, alpha]); end
end

function [phi_mod, add_pi] = phase_to_neg_half_vec(phi)
% Vectorized: map phases to (-π,0], record π flips (multiply by -1 later)
  phi = mod(phi, 2*pi);            % [0,2π)
  phi(phi > pi) = phi(phi > pi) - 2*pi;   % (-π,π]
  add_pi = phi > 0;
  phi_mod = phi;
  phi_mod(add_pi) = phi_mod(add_pi) - pi;
end

function a = phase_to_alpha(phi, w)
  TOL = 1e-12;
  if abs(phi) < TOL,        a = Inf;
  elseif abs(phi+pi) < TOL, a = 0;
  else
    a = -w * cot(phi/2);    % correct sign
    if a < 0, a = 0; end
  end
end

function [m,n] = size_local(G)
  if isnumeric(G), [m,n] = size(G);
  else, sz = size(G); m = sz(1); n = sz(2); end
end
