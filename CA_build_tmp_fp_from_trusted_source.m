function F_tmp = CA_build_tmp_fp_from_trusted_source(ap_xy, grid_xy, source_xy, Rbar_dB, n_hat, cfg)
% CA_BUILD_TMP_FP_FROM_TRUSTED_SOURCE  Propagate one trusted RSS point to grid.

    if nargin < 6 || isempty(cfg)
        cfg = CA_trusted_rss_defaults();
    elseif ~isfield(cfg, 'eps_dist') || isempty(cfg.eps_dist)
        d = CA_trusted_rss_defaults();
        cfg.eps_dist = d.eps_dist;
    end

    ap_xy = double(ap_xy);
    grid_xy = double(grid_xy);
    source_xy = double(source_xy(:)');
    Rbar_dB = double(Rbar_dB(:)');

    M = size(ap_xy, 1);
    G = size(grid_xy, 1);

    if size(ap_xy, 2) ~= 2 || size(grid_xy, 2) ~= 2 || numel(source_xy) ~= 2
        error('[CA] Invalid geometry dimensions');
    end
    if numel(Rbar_dB) ~= M
        error('[CA] Rbar_dB length must match AP count');
    end
    if ~isscalar(n_hat) || ~isfinite(n_hat) || n_hat <= 0
        error('[CA] n_hat must be a positive finite scalar');
    end

    F_tmp = nan(M, G);
    for m = 1:M
        d_mc = norm(source_xy - ap_xy(m, :));
        d_mc = max(d_mc, cfg.eps_dist);

        dx = grid_xy(:, 1) - ap_xy(m, 1);
        dy = grid_xy(:, 2) - ap_xy(m, 2);
        d_mg = sqrt(dx.^2 + dy.^2);
        d_mg = max(d_mg, cfg.eps_dist);

        F_tmp(m, :) = Rbar_dB(m) - 10 * n_hat * log10(d_mg(:)' / d_mc);
    end
end
