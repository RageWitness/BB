function SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band)
% PRECOMPUTE_CENTERED_FP_WKNN  预计算 WKNN 功率修正所需的均值和中心化指纹
%
%   SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band)
%
%   对已有的 SpatialFP_band.F_dBm (M x G)，计算：
%     mean_dBm(g)      = mean over APs of F_dBm(:,g)
%     centered_dBm(:,g) = F_dBm(:,g) - mean_dBm(g)
%
%   输入/输出：
%       SpatialFP_band - 单频带指纹结构体，需包含 F_dBm (M x G)
%       新增字段：.mean_dBm (1 x G), .centered_dBm (M x G)

    F = SpatialFP_band.F_dBm;  % M x G

    % AP 维度均值
    mu = mean(F, 1);           % 1 x G

    SpatialFP_band.mean_dBm     = mu;
    SpatialFP_band.centered_dBm = F - mu;  % 广播减法 M x G
end
