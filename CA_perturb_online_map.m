function OnlineSpatialFP = CA_perturb_online_map(OnlineSpatialFP, Config)
% CA_PERTURB_ONLINE_MAP  对在线指纹库注入空间相关扰动（模拟地图老化）
%
%   在 centered_dBm 上叠加一个 GP 采样的空间场 + per-AP 独立偏置，
%   使离线指纹库与真实信道产生可被 CA 纠正的偏差。
%
%   扰动模型：
%     δ(m,g) = A * f(g) + ε_m
%   其中:
%     f(g) ~ GP(0, SE(length_scale_m))    空间相关场（所有 AP 共享）
%     ε_m  ~ N(0, per_ap_std_dB^2)        per-AP 独立偏置
%     A    = amplitude_dB                  总幅度缩放
%
%   扰动施加到所有 band 的 centered_dBm 上。

    Config = CA_fill_defaults(Config);
    pcfg = Config.ca.perturb;

    if ~pcfg.enable
        return;
    end

    rng_state = rng;
    rng(pcfg.seed);

    grid_xy = OnlineSpatialFP.grid_xy;
    G = size(grid_xy, 1);
    M = OnlineSpatialFP.M;
    B = OnlineSpatialFP.B;
    ls = pcfg.length_scale_m;
    amp = pcfg.amplitude_dB;
    ap_std = pcfg.per_ap_std_dB;

    % --- 空间相关场 f(g) via GP ---
    D2 = zeros(G, G);
    for i = 1:G
        diff = grid_xy - grid_xy(i, :);
        D2(i, :) = sum(diff.^2, 2)';
    end
    K = exp(-0.5 * D2 / ls^2);
    K = K + 1e-4 * eye(G);
    L = chol(K, 'lower');
    f_spatial = L * randn(G, 1);
    f_spatial = f_spatial / std(f_spatial);  % normalize to unit std

    % --- per-AP 偏置 ---
    ap_bias = ap_std * randn(M, 1);

    % --- 叠加到每个 band ---
    for b = 1:B
        delta = amp * repmat(f_spatial', M, 1) + repmat(ap_bias, 1, G);
        OnlineSpatialFP.band(b).centered_dBm = ...
            OnlineSpatialFP.band(b).centered_dBm + delta;
    end

    rng(rng_state);

    OnlineSpatialFP.perturb_applied = true;
    OnlineSpatialFP.perturb_config  = pcfg;

    fprintf('[CA-Perturb] 已注入扰动: amp=%.1f dB, ls=%.0f m, ap_std=%.1f dB\n', ...
        amp, ls, ap_std);
end
