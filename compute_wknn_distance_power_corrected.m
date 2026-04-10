function [dist_vec, F_obs_centered] = compute_wknn_distance_power_corrected( ...
    F_obs, SpatialFP_band)
% COMPUTE_WKNN_DISTANCE_POWER_CORRECTED  功率修正后的 WKNN 距离计算
%
%   [dist_vec, F_obs_centered] = compute_wknn_distance_power_corrected(
%       F_obs, SpatialFP_band)
%
%   算法：
%     1. 计算观测均值：mu_obs = mean(F_obs)
%     2. 观测中心化：F_obs_centered = F_obs - mu_obs
%     3. 对每个参考点 g，计算：
%        d_g = ||F_obs_centered - centered_dBm(:,g)||_2
%
%   输入：
%       F_obs          - (M x 1) 事件平均 RSS 指纹 (dBm)
%       SpatialFP_band - 单频带指纹结构体，需含 centered_dBm (M x G)
%
%   输出：
%       dist_vec       - (1 x G) 与每个参考点的功率修正欧氏距离
%       F_obs_centered - (M x 1) 中心化后的观测指纹

    M = numel(F_obs);
    G = size(SpatialFP_band.centered_dBm, 2);

    % 1. 观测均值
    mu_obs = mean(F_obs);

    % 2. 观测中心化
    F_obs_centered = F_obs - mu_obs;  % M x 1

    % 3. 逐参考点计算距离（向量化）
    % diff = F_obs_centered - centered_dBm(:,g)，对所有 g 同时计算
    diff_mat = SpatialFP_band.centered_dBm - F_obs_centered;  % M x G
    dist_vec = sqrt(sum(diff_mat.^2, 1));                      % 1 x G
end
