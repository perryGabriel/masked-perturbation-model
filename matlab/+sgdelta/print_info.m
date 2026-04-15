function print_info(M_0, Delta_0, wstar_0, sigmax_0, info_0)
    % This function prints various computed metrics related to M, Delta, and their properties.

    % Print ||M||_inf
    fprintf('||M||_inf ≈ %.4f at ω* = %.4f rad/s\n', sigmax_0, wstar_0);

    % After building Delta with the corrected alpha/beta
    Hstar_0 = squeeze(freqresp(Delta_0, wstar_0));
    target_0 = (1/sigmax_0) * (info_0.v * conj(info_0.u.'));
    fprintf('‖Delta(jω*) - (1/σ) v u*‖₂ = %.3e  (should be ~ 0)\n', norm(Hstar_0 - target_0, 2));

    % Check the return-difference singularity
    Grw_star_0 = squeeze(freqresp(M_0, wstar_0));   % M == G_{rw}
    Ret_0 = eye(size(Grw_star_0, 1)) - Grw_star_0 * Hstar_0;
    smin_0 = min(svd(Ret_0));
    fprintf('σ_min(I - M(jω*)Δ(jω*)) = %.3e  (should be ~ 0)\n', smin_0);
end
