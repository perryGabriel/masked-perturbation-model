function [M, Grw, G] = build_M_from_ss(P, idx)
% BUILD_M_FROM_SS  Extract M = G_rw from a block-partitioned LTI.
%
% Inputs
%   P   : ss/tf/zpk LTI with inputs ordered [u; w] and outputs ordered [y; r]
%   idx : struct with fields
%           .u  (columns of inputs for u)
%           .w  (columns of inputs for w)
%           .y  (rows of outputs for y)
%           .r  (rows of outputs for r)
%
% Outputs
%   M   : G_rw (transfer from w to r)
%   Grw : same as M (for clarity)
%   G   : struct with blocks Gyu,Gyw,Gru,Grw

  if ~all(isfield(idx, {'u','w','y','r'}))
    error('idx must contain fields .u, .w, .y, .r');
  end

  Gyu = P(idx.y, idx.u);
  Gyw = P(idx.y, idx.w);
  Gru = P(idx.r, idx.u);
  Grw = P(idx.r, idx.w);

  if size(Grw,1) ~= size(Grw,2)
    warning('G_{rw} is %dx%d (non-square). The Δ loop I - G_{rw}Δ assumes square.', size(Grw,1), size(Grw,2));
  end

  M   = Grw;
  if nargout >= 3
    G = struct('Gyu',Gyu,'Gyw',Gyw,'Gru',Gru,'Grw',Grw);
  end
end
