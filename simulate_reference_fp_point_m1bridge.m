function [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge( ...
    grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
% SIMULATE_REFERENCE_FP_POINT_M1BRIDGE  在参考点生成名义 RSS 指纹
%
%   [rss_dBm, rss_lin] = simulate_reference_fp_point_m1bridge(
%       grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
%
%   在参考点 grid_pos_xy 放置一个参考单源（发射功率 ref_power_dBm），
%   调用 M1 的名义信道模型（无噪声、无慢变、无阴影随机项）生成各 AP RSS。
%
%   输入：
%       grid_pos_xy   - (1 x 2) 参考点坐标
%       band_id       - 频带编号
%       ref_power_dBm - 参考发射功率 (dBm)
%       APs           - AP 结构体
%       Bands         - 频带结构体
%       Config        - 配置（内部临时关闭慢变和阴影）
%
%   输出：
%       rss_dBm       - (M x 1) 各 AP 接收 RSS (dBm)
%       rss_lin       - (M x 1) 各 AP 接收 RSS (线性 mW)

    b = band_id;
    M = APs.num;

    % 1. 计算参考点到各 AP 距离
    dist_m = sqrt(sum((APs.pos_xy - grid_pos_xy).^2, 2));  % M x 1
    dist_m = max(dist_m, 1);  % 最小 1m 保护

    % 2. 根据频带模型计算确定性路径损耗（无阴影随机项）
    model = Bands.model{b};
    switch model
        case 'cost231wi'
            PL_dB = compute_pathloss_cost231wi(dist_m, Bands.fc_Hz(b), ...
                Config.m1.channel.cost231wi);

        case 'lognormal'
            % 确定性路径损耗：PL0 + 10*n*log10(d/d0)，不加阴影项
            d0    = Config.m1.channel.lognormal.d0_m;
            PL0_b = Config.m1.channel.lognormal.PL0_dB(b);
            n_b   = Config.m1.channel.lognormal.n_exp(b);
            dist_safe = max(dist_m, d0);
            PL_dB = PL0_b + 10 * n_b * log10(dist_safe / d0);

        otherwise
            error('[M2.5] 未知信道模型: %s', model);
    end

    % 3. 接收功率 = 发射功率 - 路径损耗（无漂移、无噪声）
    rss_dBm = ref_power_dBm - PL_dB;
    rss_lin = 10.^(rss_dBm / 10);
end
