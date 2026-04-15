function [Mmask, rcMask, idxKeep] = mask_M(M, rmask, cmask, doShrink)
% MASK_M  Zero-out selected rows/cols of M via binary masks.
%   Mmask = diag(rmask)*M*diag(cmask)   (LTI)
%   Mmask = M .* (rmask*cmask.')        (numeric)
%
% Inputs:
%   M      : ss/tf/zpk/frd (LTI) or numeric matrix
%   rmask  : |r|-by-1 (or 1-by-|r|) binary/logical vector (rows to keep=1)
%   cmask  : |w|-by-1 (or 1-by-|w|) binary/logical vector (cols to keep=1)
%   doShrink (optional, default=false): if true, also return the compressed
%            submatrix M(rmask==1, cmask==1) and the kept indices.
%
% Outputs:
%   Mmask  : masked object, same size as M (unless doShrink=true and M is numeric)
%   rcMask : outer-product mask rmask*cmask.' (numeric)
%   idxKeep: struct with fields .r (row idx kept), .c (col idx kept)

  if nargin < 4, doShrink = false; end
  rmask = logical(rmask(:));
  cmask = logical(cmask(:));
  nr = numel(rmask);  nc = numel(cmask);

  if isnumeric(M)
    % Numeric: Hadamard by outer mask
    rcMask = double(rmask) * double(cmask).';
    Mmask  = M .* rcMask;
    if doShrink
      idxKeep.r = find(rmask); idxKeep.c = find(cmask);
      Mmask = Mmask(idxKeep.r, idxKeep.c);   % compressed numeric submatrix
    else
      idxKeep = struct('r', find(rmask), 'c', find(cmask));
    end
    return;
  end

  % LTI case (ss/tf/zpk/frd): pre/post multiply by diagonal masks
  Dl = diag(double(rmask));   % size nr x nr
  Dr = diag(double(cmask));   % size nc x nc
  rcMask = double(rmask) * double(cmask).';

  % Safety: check I/O sizes
  [p,q] = size(M);
  if p ~= nr || q ~= nc
    error('Mask sizes do not match M: M is %dx%d, rmask=%d, cmask=%d.', p,q,nr,nc);
  end

  Mmask = Dl * M * Dr;

  % Optional: shrink by selecting channels (keeps LTI typing intact)
  if doShrink
    idxKeep.r = find(rmask);
    idxKeep.c = find(cmask);
    % Channel selection (same as taking submatrix of I/O)
    Mmask = Mmask(idxKeep.r, idxKeep.c);
  else
    idxKeep = struct('r', find(rmask), 'c', find(cmask));
  end
end
