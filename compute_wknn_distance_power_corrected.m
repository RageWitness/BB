function [dist_vec, F_obs_centered] = compute_wknn_distance_power_corrected( ...
    F_obs, SpatialFP_band)
% COMPUTE_WKNN_DISTANCE_POWER_CORRECTED  功率修正后的 WKNN 距离计算
%
%   [dist_vec, F_obs_centered] = compute_wknn_distance_power_corrected(
%       F_obs, SpatialFP_band)
%
%   算法（双项距离）：
%     d_g^2 = lambda_c * ||F_obs_c - F_db_c(:,g)||^2
%           + lambda_mu * (mean_db(g) - mu_obs - alpha0)^2
%
%     第一项：形状距离（中心化后欧氏距离）
%     第二项：均值层级惩罚（保留绝对功率信息）
%
%   输入：
%       F_obs          - (M x 1) 事件平均 RSS 指纹 (dBm)
%       SpatialFP_band - 单频带指纹结构体，需含：
%                          centered_dBm (M x G) 中心化指纹
%                          mean_dBm     (1 x G) 各网格点均值
%
%   输出：
%       dist_vec       - (1 x G) 与每个参考点的距离
%       F_obs_centered - (M x 1) 中心化后的观测指纹

    % ---- 可调参数 ----
    lambda_c  = 1.0;   % 形状距离权重
    lambda_mu = 1.0;   % 均值层级惩罚权重
    alpha0    = 0;     % 库参考功率与目标发射功率之差 (dB)
    % ------------------

    % 1. 观测均值 & 中心化
    mu_obs = mean(F_obs);
    F_obs_centered = F_obs - mu_obs;  % M x 1

    % 2. 形状距离（中心化欧氏）
    diff_shape = SpatialFP_band.centered_dBm - F_obs_centered;  % M x G
    d_shape_sq = sum(diff_shape.^2, 1);                          % 1 x G

    % 3. 均值层级惩罚
    d_mu_sq = (SpatialFP_band.mean_dBm - mu_obs - alpha0).^2;   % 1 x G

    % 4. 合成距离
    dist_vec = sqrt(lambda_c * d_shape_sq + lambda_mu * d_mu_sq);
end
