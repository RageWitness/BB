function SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band)
% PRECOMPUTE_CENTERED_FP_WKNN  预计算 WKNN 所需的均值和中心化指纹
%
%   SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band)
%
%   dBm 域：
%     mean_dBm(g)       = mean over APs of F_dBm(:,g)
%     centered_dBm(:,g) = F_dBm(:,g) - mean_dBm(g)    （减法去均值）
%
%   线性域（用于线性域匹配）：
%     mean_lin(g)       = mean over APs of F_lin(:,g)
%     centered_lin(:,g) = F_lin(:,g) / mean_lin(g)     （除法归一化）
%
%   输入/输出：
%       SpatialFP_band - 单频带指纹结构体，需包含 F_dBm (M x G), F_lin (M x G)
%       新增字段：
%         .mean_dBm (1 x G), .centered_dBm (M x G)
%         .mean_lin (1 x G), .centered_lin (M x G)

    % --- dBm 域 ---
    F_dBm = SpatialFP_band.F_dBm;       % M x G
    mu_dBm = mean(F_dBm, 1);            % 1 x G
    SpatialFP_band.mean_dBm     = mu_dBm;
    SpatialFP_band.centered_dBm = F_dBm - mu_dBm;

    % --- 线性域 ---
    F_lin = SpatialFP_band.F_lin;        % M x G
    mu_lin = mean(F_lin, 1);             % 1 x G
    mu_lin_safe = max(mu_lin, 1e-30);    % 防除零
    SpatialFP_band.mean_lin     = mu_lin;
    SpatialFP_band.centered_lin = F_lin ./ mu_lin_safe;  % 除法归一化
end
