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
        ref_power = 20 * ones(1, B);
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

    % --- 判断是否需要几何缓存（MWM 模型才需要）---
    need_geometry = false;
    for b = 1:B
        if strcmp(Bands.model{b}, 'mwm')
            need_geometry = true;
            break;
        end
    end

    GeometryCache = [];
    if need_geometry
        buildings = Config.m1.channel.buildings;
        GeometryCache = precompute_geometry_cache(GridValid.xy, APs.pos_xy, buildings);
    end

    N_mc = 50;
    if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
        N_mc = Config.m25.n_mc_fp;
    end

    % --- 噪声功率（每频带不同带宽） ---
    N_dBm_per_band = zeros(1, B);
    if isfield(Config, 'm1') && isfield(Config.m1, 'noise')
        noise_cfg = Config.m1.noise;
        if isfield(noise_cfg, 'mode') && strcmp(noise_cfg.mode, 'power_dBm') ...
                && isfield(noise_cfg, 'noise_power_dBm')
            for b = 1:B
                N_dBm_per_band(b) = noise_cfg.noise_power_dBm(b);
            end
        else
            n0_arr = -174 * ones(1, B);
            if isfield(noise_cfg, 'n0_dBmHz')
                n0_arr = noise_cfg.n0_dBmHz;
                if isscalar(n0_arr), n0_arr = n0_arr * ones(1, B); end
            end
            for b = 1:B
                N_dBm_per_band(b) = n0_arr(b) + 10 * log10(Bands.bw_Hz(b));
            end
        end
    else
        for b = 1:B
            N_dBm_per_band(b) = -174 + 10 * log10(Bands.bw_Hz(b));
        end
    end
    N_lin_per_band = 10.^(N_dBm_per_band / 10);

    for b = 1:B
        F_dBm   = zeros(M, G);
        F_lin   = zeros(M, G);
        Std_dBm = zeros(M, G);
        Std_lin = zeros(M, G);

        N_lin_b = N_lin_per_band(b);

        if strcmp(Bands.model{b}, 'lognormal')
            % --- 对数正态阴影衰落模型 ---
            ln = Config.m1.channel.lognormal;
            d0_ln    = ln.d0;
            n_exp    = ln.n;
            sigma_b  = ln.sigma;
            if isfield(ln, 'PL0')
                PL0_b = ln.PL0;
            else
                PL0_b = 20 * log10(Bands.fc_Hz(b) / 1e6) - 28;
            end

            for g = 1:G
                dist_g = sqrt(sum((APs.pos_xy - GridValid.xy(g,:)).^2, 2));
                dist_g = max(dist_g, d0_ln);
                PL_det = PL0_b + 10 * n_exp * log10(dist_g / d0_ln);

                xi = sigma_b * randn(M, N_mc);
                rss_dBm_samp = ref_power(b) - PL_det + xi;
                rss_lin_samp = 10.^(rss_dBm_samp / 10);

                % Add the mean receiver noise power floor in mW.
                % N_lin_b is already receiver noise power.
                rss_lin_samp = max(rss_lin_samp + N_lin_b, 1e-20);

                rss_lin_mean = mean(rss_lin_samp, 2);
                F_dBm(:, g)   = 10 * log10(max(rss_lin_mean, 1e-20));
                F_lin(:, g)   = rss_lin_mean;
                Std_dBm(:, g) = std(10 * log10(max(rss_lin_samp, 1e-20)), 0, 2);
                Std_lin(:, g) = std(rss_lin_samp, 0, 2);
            end
        else
            % --- MWM 模型 ---
            p_band  = get_channel_params(Bands.fc_Hz(b));
            sigma_b = p_band.sigma;

            for g = 1:G
                PL_det = compute_pathloss_from_cache(g, Bands.fc_Hz(b), GeometryCache);
                xi = sigma_b * randn(M, N_mc);
                rss_dBm_samp = ref_power(b) - PL_det + xi;
                rss_lin_samp = 10.^(rss_dBm_samp / 10);

                % Add the mean receiver noise power floor in mW.
                % N_lin_b is already receiver noise power.
                rss_lin_samp = max(rss_lin_samp + N_lin_b, 1e-20);

                rss_lin_mean = mean(rss_lin_samp, 2);
                F_dBm(:, g)   = 10 * log10(max(rss_lin_mean, 1e-20));
                F_lin(:, g)   = rss_lin_mean;
                Std_dBm(:, g) = std(10 * log10(max(rss_lin_samp, 1e-20)), 0, 2);
                Std_lin(:, g) = std(rss_lin_samp, 0, 2);
            end
        end

        band_fp.F_dBm   = F_dBm;
        band_fp.F_lin   = F_lin;
        band_fp.std_dBm = Std_dBm;
        band_fp.std_lin = Std_lin;
        band_fp.fc_Hz   = Bands.fc_Hz(b);
        band_fp.bw_Hz   = Bands.bw_Hz(b);
        band_fp.model   = Bands.model{b};
        band_fp.noise_floor_dBm = N_dBm_per_band(b);
        band_fp.noise_model = 'mean_power_floor';
        band_fp.reference_power_dBm_used = ref_power(b);

        % 预计算 RF_raw / RF_minmax / F_shape_l1 / centered_dBm
        band_fp = precompute_centered_fp_wknn(band_fp, fp_mode);

        SpatialFP.band(b) = band_fp;

        fprintf('  Band %d (%s): N_floor=%.1f dBm, F_dBm [%.1f, %.1f], F_lin [%.2e, %.2e]', ...
            b, Bands.name{b}, N_dBm_per_band(b), min(F_dBm(:)), max(F_dBm(:)), ...
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
