function [PL_dB, sigma_dB] = compute_pathloss_lognormal(src_xy, ap_xy_all, freq_hz, Config)
% COMPUTE_PATHLOSS_LOGNORMAL  对数正态阴影衰落路径损耗模型
%
%   PL(d) = PL_0 + 10*n*log10(d/d0) + X_sigma
%
%   其中 X_sigma ~ N(0, sigma^2)，每次调用独立采样。
%
%   输入：
%     src_xy     - (1 x 2) 源坐标
%     ap_xy_all  - (M x 2) AP 坐标
%     freq_hz    - 载频 Hz
%     Config     - 配置（读取 Config.m1.channel.lognormal）
%
%   输出：
%     PL_dB    - (M x 1) 含阴影衰落的路径损耗
%     sigma_dB - 标量，阴影衰落标准差

    M = size(ap_xy_all, 1);

    ln = Config.m1.channel.lognormal;
    d0    = ln.d0;
    n_exp = ln.n;
    sigma = ln.sigma;

    if isfield(ln, 'PL0')
        PL0 = ln.PL0;
    else
        f_MHz = freq_hz / 1e6;
        PL0 = 20 * log10(f_MHz) - 28;
    end

    dist = sqrt(sum((ap_xy_all - src_xy).^2, 2));
    dist = max(dist, d0);

    PL_det = PL0 + 10 * n_exp * log10(dist / d0);

    X_sigma = sigma * randn(M, 1);

    PL_dB = PL_det + X_sigma;
    sigma_dB = sigma;
end
