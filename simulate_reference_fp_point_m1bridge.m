function [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin, rss_lin_samples] = simulate_reference_fp_point_m1bridge( ...
    grid_pos_xy, band_id, ref_power_dBm, APs, Bands, Config)
% SIMULATE_REFERENCE_FP_POINT_M1BRIDGE  在参考点生成名义 RSS 指纹
%
%   统一 MWM 模型：确定性路径损耗 + N(0, sigma^2) 残差 MC 采样
%   建库时取线性域均值作为指纹

    b = band_id;
    M = APs.num;

    buildings = Config.m1.channel.buildings;

    % 确定性 MWM 路径损耗（按 AP 逐个）
    [PL_det, sigma_b] = compute_pathloss_mwm(grid_pos_xy, APs.pos_xy, Bands.fc_Hz(b), buildings);

    % MC 样本数
    N_mc = 50;
    if isfield(Config, 'm25') && isfield(Config.m25, 'n_mc_fp')
        N_mc = Config.m25.n_mc_fp;
    end

    rss_dBm_samples = zeros(M, N_mc);
    rss_lin_samples = zeros(M, N_mc);
    for mc = 1:N_mc
        xi = sigma_b * randn(M, 1);
        rss_mc_dBm = ref_power_dBm - PL_det + xi;
        rss_dBm_samples(:, mc) = rss_mc_dBm;
        rss_lin_samples(:, mc) = 10.^(rss_mc_dBm / 10);
    end

    rss_lin = mean(rss_lin_samples, 2);
    rss_dBm = 10 * log10(max(rss_lin, 1e-20));
    rss_std_dBm = std(rss_dBm_samples, 0, 2);
    rss_std_lin = std(rss_lin_samples, 0, 2);
end
