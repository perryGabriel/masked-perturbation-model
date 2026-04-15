function [M_i, Delta_i, wstar_i, sigmax_i, info_i] = process_masks(M, rmask, wmask, wgrid)
    % This function processes the masks and returns M_i, Delta_i, and info_i.

    % Apply the mask to M
    [M_i, ~] = mask_M(M, rmask, wmask, false);

    % Create Delta and info using the modified M_i
    [Delta_i, wstar_i, sigmax_i, info_i] = construct_delta_allpass(M_i, wgrid, struct('verbose', false));
end
