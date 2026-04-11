function [rss_dBm, rss_lin, rss_std_dBm] = simulate_reference_fp_point_m1bridge( ...
    grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
% SIMULATE_REFERENCE_FP_POINT_M1BRIDGE  在参考点生成名义 RSS 指纹
%
%   [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge(...)
%   [rss_dBm, rss_lin, rss_std_dBm] = simulate_reference_fp_point_m1bridge(...)
%
%   对 lognormal 频带：进行 N_mc 次 Monte Carlo 采样（含 i.i.d. 阴影），
%   取线性域均值作为指纹，同时输出 dBm 域标准差。
%   对 cost231wi 频带：路径损耗确定性，std = 0。
%
%   输入：
%       grid_pos_xy   - (1 x 2) 参考点坐标
%       band_id       - 频带编号
%       ref_power_dBm - 参考发射功率 (dBm)
%       APs           - AP 结构体
%       Bands         - 频带结构体
%       Config        - 配置（含 Config.m25.n_mc_fp 采样次数，默认 50）
%
%   输出：
%       rss_dBm     - (M x 1) 各 AP 接收 RSS 均值 (dBm)
%       rss_lin     - (M x 1) 各 AP 接收 RSS 均值 (线性 mW)
%       rss_std_dBm - (M x 1) 各 AP RSS 标准差 (dB)，cost231wi 为 0

    b = band_id;
    M = APs.num;

    % 1. 计算参考点到各 AP 距离
    dist_m = sqrt(sum((APs.pos_xy - grid_pos_xy).^2, 2));  % M x 1
    dist_m = max(dist_m, 1);  % 最小 1m 保护

    % 2. 根据频带模型计算路径损耗
    model = Bands.model{b};
    switch model
        case 'cost231wi'
            % 确定性模型
            PL_dB = compute_pathloss_cost231wi(dist_m, Bands.fc_Hz(b), ...
                Config.m1.channel.cost231wi);
            rss_dBm     = ref_power_dBm - PL_dB;
            rss_lin     = 10.^(rss_dBm / 10);
            rss_std_dBm = zeros(M, 1);

        case 'lognormal'
            % Monte Carlo 采样
            N_mc = 50;
            if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
                N_mc = Config.m25.n_mc_fp;
            end

            d0      = Config.m1.channel.lognormal.d0_m;
            PL0_b   = Config.m1.channel.lognormal.PL0_dB(b);
            n_b     = Config.m1.channel.lognormal.n_exp(b);
            sigma_b = Config.m1.channel.lognormal.sigma_dB(b);

            dist_safe = max(dist_m, d0);
            PL_det = PL0_b + 10 * n_b * log10(dist_safe / d0);  % M x 1

            % 收集每次采样的 dBm 值
            rss_dBm_samples = zeros(M, N_mc);
            rss_lin_sum     = zeros(M, 1);
            for mc = 1:N_mc
                X_shadow = sigma_b * randn(M, 1);
                PL_mc    = PL_det + X_shadow;
                rss_mc_dBm = ref_power_dBm - PL_mc;
                rss_mc_lin = 10.^(rss_mc_dBm / 10);
                rss_dBm_samples(:, mc) = rss_mc_dBm;
                rss_lin_sum = rss_lin_sum + rss_mc_lin;
            end

            % 线性域均值 → dBm
            rss_lin = rss_lin_sum / N_mc;
            rss_dBm = 10 * log10(max(rss_lin, 1e-20));

            % dBm 域标准差
            rss_std_dBm = std(rss_dBm_samples, 0, 2);  % M x 1

        otherwise
            error('[M2.5] 未知信道模型: %s', model);
    end
end
