function SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config)
% BUILD_SPATIAL_FP_SINGLE_SOURCE  构建每个频带的单源空间指纹库
%
%   仅构建：F_dBm / F_lin / RF_raw / RF_minmax / F_shape_l1 / centered_dBm
%   （概率 shape / hard-negative 已移除）

    B = Bands.B;
    M = APs.num;
    G = GridValid.Nvalid;

    % 默认参考功率
    if isfield(Config, 'm25') && isfield(Config.m25, 'ref_power_dBm')
        ref_power = Config.m25.ref_power_dBm;
    else
        ref_power = 100 * ones(1, B);
    end

    % 指纹模式
    fp_mode = 'both';
    if isfield(Config, 'm25') && isfield(Config.m25, 'fingerprint_mode')
        fm = Config.m25.fingerprint_mode;
        if strcmp(fm, 'legacy_only'), fp_mode = 'legacy'; end
    end
    if isfield(Config, 'm25') && isfield(Config.m25, 'keep_legacy_shape') ...
            && ~logical(Config.m25.keep_legacy_shape)
        fp_mode = 'rf_minmax';
    end

    SpatialFP.B = B;
    SpatialFP.M = M;
    SpatialFP.G = G;
    SpatialFP.grid_xy       = GridValid.xy;
    SpatialFP.ref_power_dBm = ref_power;
    SpatialFP.reference_power_dBm_used = ref_power;
    SpatialFP.shape_library_mode = 'deterministic_l1_from_mean_lin';
    SpatialFP.fingerprint_mode = fp_mode;

    fprintf('[M2.5] 构建 SpatialFP: M=%d AP, G=%d 网格点, B=%d 频带\n', M, G, B);

    % --- 几何缓存（与频率无关，一次构建供 B 频带复用）---
    buildings = Config.m1.channel.buildings;
    GeometryCache = precompute_geometry_cache(GridValid.xy, APs.pos_xy, buildings);

    N_mc = 50;
    if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
        N_mc = Config.m25.n_mc_fp;
    end

    for b = 1:B
        F_dBm   = zeros(M, G);
        F_lin   = zeros(M, G);
        Std_dBm = zeros(M, G);
        Std_lin = zeros(M, G);

        p_band  = get_channel_params(Bands.fc_Hz(b));
        sigma_b = p_band.sigma;

        for g = 1:G
            PL_det = compute_pathloss_from_cache(g, Bands.fc_Hz(b), GeometryCache);
            xi = sigma_b * randn(M, N_mc);
            rss_dBm_samp = ref_power(b) - PL_det + xi;
            rss_lin_samp = 10.^(rss_dBm_samp / 10);

            rss_lin_mean = mean(rss_lin_samp, 2);
            F_dBm(:, g)   = 10 * log10(max(rss_lin_mean, 1e-20));
            F_lin(:, g)   = rss_lin_mean;
            Std_dBm(:, g) = std(rss_dBm_samp, 0, 2);
            Std_lin(:, g) = std(rss_lin_samp, 0, 2);
        end

        band_fp.F_dBm   = F_dBm;
        band_fp.F_lin   = F_lin;
        band_fp.std_dBm = Std_dBm;
        band_fp.std_lin = Std_lin;
        band_fp.fc_Hz   = Bands.fc_Hz(b);
        band_fp.bw_Hz   = Bands.bw_Hz(b);
        band_fp.model   = Bands.model{b};
        band_fp.reference_power_dBm_used = ref_power(b);

        % 预计算 RF_raw / RF_minmax / F_shape_l1 / centered_dBm
        band_fp = precompute_centered_fp_wknn(band_fp, fp_mode);

        SpatialFP.band(b) = band_fp;

        fprintf('  Band %d (%s): F_dBm [%.1f, %.1f], F_lin [%.2e, %.2e]', ...
            b, Bands.name{b}, min(F_dBm(:)), max(F_dBm(:)), ...
            min(F_lin(:)), max(F_lin(:)));
        if isfield(band_fp, 'RF_minmax')
            fprintf(', RF_minmax [%.3f, %.3f]', ...
                min(band_fp.RF_minmax(:)), max(band_fp.RF_minmax(:)));
        end
        if isfield(band_fp, 'F_shape_l1')
            fprintf(', shape_l1 [%.4f, %.4f]', ...
                min(band_fp.F_shape_l1(:)), max(band_fp.F_shape_l1(:)));
        end
        fprintf('\n');
    end

    fprintf('[M2.5] SpatialFP 构建完成\n');
end
