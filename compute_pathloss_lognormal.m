function PL_dB = compute_pathloss_lognormal(dist_m, band_id, ChannelState, Config)
% COMPUTE_PATHLOSS_LOGNORMAL  对数正态阴影路径损耗模型
%
%   PL_dB = compute_pathloss_lognormal(dist_m, band_id, ChannelState, Config)
%
%   公式：
%     PL(m,b) = PL0(b) + 10*n(b)*log10(d(m)/d0) + X(m,b)
%
%   其中 X(m,b) 为阴影项：
%     若 enableSlowDrift=true  → 使用 ChannelState.shadow_state_dB(:,b)
%     若 enableSlowDrift=false → 当场生成 i.i.d. N(0, sigma_b^2)
%
%   输入：
%       dist_m       - (M x 1) 源到各 AP 距离 (m)
%       band_id      - 频带编号
%       ChannelState - 信道状态
%       Config       - 配置
%
%   输出：
%       PL_dB        - (M x 1) 路径损耗 (dB)

    b = band_id;
    M = numel(dist_m);

    d0       = Config.m1.channel.lognormal.d0_m;
    PL0_b    = Config.m1.channel.lognormal.PL0_dB(b);
    n_b      = Config.m1.channel.lognormal.n_exp(b);
    sigma_b  = Config.m1.channel.lognormal.sigma_dB(b);

    % 距离下限保护（避免 log10(0)）
    dist_m = max(dist_m, d0);

    % 确定性路径损耗
    PL_det = PL0_b + 10 * n_b * log10(dist_m / d0);

    % 阴影项
    if Config.m1.channel.enableSlowDrift
        % 使用已更新的 AR(1) 阴影状态
        X_shadow = ChannelState.shadow_state_dB(:, b);
    else
        % i.i.d. 阴影
        X_shadow = sigma_b * randn(M, 1);
    end

    PL_dB = PL_det + X_shadow;
end
