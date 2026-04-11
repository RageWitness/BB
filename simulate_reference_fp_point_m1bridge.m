function [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge( ...
    grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
% SIMULATE_REFERENCE_FP_POINT_M1BRIDGE  在参考点生成名义 RSS 指纹
%
%   [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge(
%       grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
%
%   在参考点 grid_pos_xy 放置一个参考单源（发射功率 ref_power_dBm），
%   调用 M1 的名义信道模型生成各 AP RSS。
%
%   对 lognormal 频带：进行 N_mc 次 Monte Carlo 采样（含 i.i.d. 阴影），
%   取均值作为指纹，使库指纹包含阴影效应的统计平均。
%   对 cost231wi 频带：路径损耗确定性，无需采样。
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
%       rss_dBm       - (M x 1) 各 AP 接收 RSS (dBm)
%       rss_lin       - (M x 1) 各 AP 接收 RSS (线性 mW)

    b = band_id;
    M = APs.num;

    % 1. 计算参考点到各 AP 距离
    dist_m = sqrt(sum((APs.pos_xy - grid_pos_xy).^2, 2));  % M x 1
    dist_m = max(dist_m, 1);  % 最小 1m 保护

    % 2. 根据频带模型计算路径损耗
    model = Bands.model{b};
    switch model
        case 'cost231wi'
            % 确定性模型，无需 Monte Carlo
            PL_dB = compute_pathloss_cost231wi(dist_m, Bands.fc_Hz(b), ...
                Config.m1.channel.cost231wi);
            rss_dBm = ref_power_dBm - PL_dB;

        case 'lognormal'
            % Monte Carlo 采样：含 i.i.d. 阴影
            N_mc = 50;  % 默认采样次数
            if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
                N_mc = Config.m25.n_mc_fp;
            end

            d0      = Config.m1.channel.lognormal.d0_m;
            PL0_b   = Config.m1.channel.lognormal.PL0_dB(b);
            n_b     = Config.m1.channel.lognormal.n_exp(b);
            sigma_b = Config.m1.channel.lognormal.sigma_dB(b);

            dist_safe = max(dist_m, d0);
            PL_det = PL0_b + 10 * n_b * log10(dist_safe / d0);  % M x 1

            % 在线性域取均值（对数正态的均值在线性域更准确）
            rss_lin_sum = zeros(M, 1);
            for mc = 1:N_mc
                X_shadow = sigma_b * randn(M, 1);            % dB 域阴影
                PL_mc    = PL_det + X_shadow;                 % M x 1
                rss_mc   = 10.^((ref_power_dBm - PL_mc) / 10);  % 线性 mW
                rss_lin_sum = rss_lin_sum + rss_mc;
            end
            rss_lin_avg = rss_lin_sum / N_mc;
            rss_dBm = 10 * log10(rss_lin_avg);

        otherwise
            error('[M2.5] 未知信道模型: %s', model);
    end

    rss_lin = 10.^(rss_dBm / 10);
end
