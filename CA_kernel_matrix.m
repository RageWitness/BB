function K = CA_kernel_matrix(X1, X2, theta, band_id, ap_id, KernelCtx, Config)
% CA_KERNEL_MATRIX  CA kernel with Matern/SE and propagation-aware distance.

    Config = CA_fill_defaults(Config);
    gpr_cfg = Config.ca.gpr;

    if nargin < 6
        KernelCtx = [];
    end

    ell_x = exp(theta(1));
    ell_y = exp(theta(2));
    sf    = exp(theta(3));

    X1 = ensure_2d(X1);
    X2 = ensure_2d(X2);

    dx = (X1(:,1) - X2(:,1)').^2 / max(ell_x^2, eps);
    dy = (X1(:,2) - X2(:,2)').^2 / max(ell_y^2, eps);
    d2 = dx + dy;

    mode = lower(gpr_cfg.distance_mode);
    if ~strcmp(mode, 'geo') && ~isempty(KernelCtx)
        weights = gpr_cfg.dist_weight;
        d2 = weights.alpha_geo * d2;

        need_fp = contains(mode, 'fp');
        need_pl = contains(mode, 'pl');
        if need_fp || need_pl
            Z1 = interpolate_centered_fp(X1, band_id, KernelCtx, Config);
            Z2 = interpolate_centered_fp(X2, band_id, KernelCtx, Config);
        end

        if need_fp
            fp2 = pairwise_sqdist(Z1, Z2);
            s_fp = max(KernelCtx.band(band_id).s_fp, gpr_cfg.dist_norm.eps);
            d2 = d2 + weights.beta_fp * fp2 / (s_fp^2);
        end

        if need_pl
            z1 = Z1(:, ap_id);
            z2 = Z2(:, ap_id);
            pl2 = (z1 - z2').^2;
            s_pl = max(KernelCtx.band(band_id).s_pl(ap_id), gpr_cfg.dist_norm.eps);
            d2 = d2 + weights.gamma_pl * pl2 / (s_pl^2);
        end
    end

    d2 = max(d2, 0);
    d = sqrt(d2);

    switch lower(gpr_cfg.kernel_type)
        case 'se'
            K = sf^2 * exp(-0.5 * d2);
        case 'matern32'
            q = sqrt(3) * d;
            K = sf^2 * (1 + q) .* exp(-q);
        case 'matern52'
            q = sqrt(5) * d;
            K = sf^2 * (1 + q + (q.^2)/3) .* exp(-q);
        otherwise
            error('[CA] Unknown kernel_type=%s', gpr_cfg.kernel_type);
    end
end


function X = ensure_2d(X)
    if isempty(X)
        X = zeros(0, 2);
    elseif size(X, 2) ~= 2 && size(X, 1) == 2
        X = X';
    end
end


function D2 = pairwise_sqdist(A, B)
    aa = sum(A.^2, 2);
    bb = sum(B.^2, 2)';
    D2 = aa + bb - 2 * (A * B');
    D2 = max(D2, 0);
end


function Zx = interpolate_centered_fp(X, band_id, KernelCtx, Config)
    Zmap = KernelCtx.band(band_id).centered_dBm;  % M x G
    grid_xy = KernelCtx.grid_xy;
    interp_cfg = Config.ca.interp;
    eps_v = interp_cfg.eps;
    K_i = interp_cfg.k;
    p_exp = interp_cfg.power;

    N = size(X, 1);
    M = size(Zmap, 1);
    Zx = zeros(N, M);

    for i = 1:N
        key = sprintf('%.9g_%.9g', X(i, 1), X(i, 2));
        if isKey(KernelCtx.grid_index_map, key)
            idx = KernelCtx.grid_index_map(key);
            Zx(i, :) = Zmap(:, idx)';
            continue;
        end

        dists = sqrt(sum((grid_xy - X(i, :)).^2, 2));
        [dmin, idx_min] = min(dists);
        if dmin < eps_v
            Zx(i, :) = Zmap(:, idx_min)';
            continue;
        end

        [~, sort_idx] = sort(dists);
        nn_idx = sort_idx(1:min(K_i, numel(sort_idx)));
        d_nn = dists(nn_idx);
        w = 1 ./ (d_nn .^ p_exp + eps_v);
        w = w / sum(w);
        Zx(i, :) = (Zmap(:, nn_idx) * w)';
    end
end
