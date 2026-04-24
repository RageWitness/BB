function KernelCtx = CA_prepare_kernel_context(OnlineSpatialFP, Config)
% CA_PREPARE_KERNEL_CONTEXT  Prepare scales for propagation-aware CA kernels.
%
%   The context stores centered maps and robust distance scales. Pairwise
%   scales are estimated from sampled grid pairs; no G-by-G matrices are built.

    Config = CA_fill_defaults(Config);
    gpr_cfg = Config.ca.gpr;
    n_pair = gpr_cfg.dist_norm.n_pair_sample;

    KernelCtx = struct();
    KernelCtx.grid_xy = OnlineSpatialFP.grid_xy;
    KernelCtx.B = OnlineSpatialFP.B;
    KernelCtx.M = OnlineSpatialFP.M;
    KernelCtx.G = OnlineSpatialFP.G;
    KernelCtx.kernel_type = gpr_cfg.kernel_type;
    KernelCtx.distance_mode = gpr_cfg.distance_mode;
    KernelCtx.alpha_geo = gpr_cfg.dist_weight.alpha_geo;
    KernelCtx.beta_fp = gpr_cfg.dist_weight.beta_fp;
    KernelCtx.gamma_pl = gpr_cfg.dist_weight.gamma_pl;

    keys = cell(OnlineSpatialFP.G, 1);
    for g = 1:OnlineSpatialFP.G
        keys{g} = make_grid_key(OnlineSpatialFP.grid_xy(g, :));
    end
    KernelCtx.grid_index_map = containers.Map(keys, num2cell((1:OnlineSpatialFP.G)'));

    for b = 1:OnlineSpatialFP.B
        if ~isfield(OnlineSpatialFP.band(b), 'centered_dBm')
            error('[CA] Kernel context requires SpatialFP.band(%d).centered_dBm', b);
        end

        Z = OnlineSpatialFP.band(b).centered_dBm;
        KernelCtx.band(b).centered_dBm = Z;

        G = size(Z, 2);
        if G < 2
            KernelCtx.band(b).s_fp = 1;
            KernelCtx.band(b).s_pl = ones(1, size(Z, 1));
            continue;
        end

        np = min(max(1, n_pair), G * max(G - 1, 1));
        idx1 = randi(G, np, 1);
        idx2 = randi(G, np, 1);
        same = idx1 == idx2;
        idx2(same) = mod(idx2(same), G) + 1;

        d_fp = sqrt(sum((Z(:, idx1) - Z(:, idx2)).^2, 1));
        KernelCtx.band(b).s_fp = robust_scale(d_fp);

        M = size(Z, 1);
        s_pl = ones(1, M);
        for m = 1:M
            d_pl = abs(Z(m, idx1) - Z(m, idx2));
            s_pl(m) = robust_scale(d_pl);
        end
        KernelCtx.band(b).s_pl = s_pl;
    end
end


function s = robust_scale(v)
    v = v(:);
    v = v(isfinite(v) & v > 0);
    if isempty(v)
        s = 1;
    else
        s = median(v);
        if ~isfinite(s) || s <= 0
            s = 1;
        end
    end
end


function key = make_grid_key(xy)
    key = sprintf('%.9g_%.9g', xy(1), xy(2));
end
