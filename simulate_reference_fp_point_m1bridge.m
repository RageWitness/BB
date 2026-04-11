function [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin] = simulate_reference_fp_point_m1bridge( ...
    grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
% SIMULATE_REFERENCE_FP_POINT_M1BRIDGE  在参考点生成名义 RSS 指纹
%
%   [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge(...)
%   [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin] = simulate_reference_fp_point_m1bridge(...)
%
%   对 lognormal 频带：N_mc 次 Monte Carlo（含 i.i.d. 阴影），
%   取线性域均值作为指纹，同时输出 dBm 域和线性域标准差。
%   对 cost231wi 频带：确定性路径损耗，std = 0。
%
%   输入：
%       grid_pos_xy   - (1 x 2) 参考点坐标
%       band_id       - 频带编号
%       ref_power_dBm - 参考发射功率 (dBm)
%       APs, Bands, Config
%
%   输出：
%       rss_dBm     - (M x 1) RSS 均值 (dBm)
%       rss_lin     - (M x 1) RSS 均值 (线性 mW)
%       rss_std_dBm - (M x 1) RSS 标准差 (dB)
%       rss_std_lin - (M x 1) RSS 标准差 (线性 mW)

    b = band_id;
    M = APs.num;

    % 1. 距离
    dist_m = sqrt(sum((APs.pos_xy - grid_pos_xy).^2, 2));
    dist_m = max(dist_m, 1);

    % 2. 路径损耗
    model = Bands.model{b};
    switch model
        case 'cost231wi'
            PL_dB = compute_pathloss_cost231wi(dist_m, Bands.fc_Hz(b), ...
                Config.m1.channel.cost231wi);
            rss_dBm     = ref_power_dBm - PL_dB;
            rss_lin     = 10.^(rss_dBm / 10);
            rss_std_dBm = zeros(M, 1);
            rss_std_lin = zeros(M, 1);

        case 'lognormal'
            N_mc = 50;
            if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
                N_mc = Config.m25.n_mc_fp;
            end

            d0      = Config.m1.channel.lognormal.d0_m;
            PL0_b   = Config.m1.channel.lognormal.PL0_dB(b);
            n_b     = Config.m1.channel.lognormal.n_exp(b);
            sigma_b = Config.m1.channel.lognormal.sigma_dB(b);

            dist_safe = max(dist_m, d0);
            PL_det = PL0_b + 10 * n_b * log10(dist_safe / d0);

            rss_dBm_samples = zeros(M, N_mc);
            rss_lin_samples = zeros(M, N_mc);
            for mc = 1:N_mc
                X_shadow = sigma_b * randn(M, 1);
                PL_mc    = PL_det + X_shadow;
                rss_mc_dBm = ref_power_dBm - PL_mc;
                rss_mc_lin = 10.^(rss_mc_dBm / 10);
                rss_dBm_samples(:, mc) = rss_mc_dBm;
                rss_lin_samples(:, mc) = rss_mc_lin;
            end

            % 线性域均值
            rss_lin = mean(rss_lin_samples, 2);
            rss_dBm = 10 * log10(max(rss_lin, 1e-20));

            % 两域标准差
            rss_std_dBm = std(rss_dBm_samples, 0, 2);
            rss_std_lin = std(rss_lin_samples, 0, 2);

        otherwise
            error('[M2.5] 未知信道模型: %s', model);
    end
end
