function [D_vec, d_shape_vec, resid_vec, q_opt_vec] = compute_wknn_distance_shape_scale( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_SHAPE_SCALE  基于形状+缩放的 WKNN 距离
%
%   [D_vec, d_shape_vec, resid_vec, q_opt_vec] = compute_wknn_distance_shape_scale(
%       F_obs_lin, F_obs_shape_l1, SpatialFP_band)
%   [D_vec, d_shape_vec, resid_vec, q_opt_vec] = compute_wknn_distance_shape_scale(
%       F_obs_lin, F_obs_shape_l1, SpatialFP_band, match_cfg)
%
%   观测模型：y ≈ q * phi_g + n
%     y      — 线性功率观测向量 (M x 1)
%     phi_g  — 参考点 g 的线性功率指纹 (M x 1)
%     q      — 未知功率缩放因子
%
%   算法：
%     1. 形状距离 (L2)：
%           d_shape(g) = ||s_obs - s_g||_2
%        其中 s = F_lin / (||F_lin||_1 + eps) 为 L1 归一化形状
%
%     2. 最优缩放系数（最小二乘解析解）：
%           q_g* = max(0, phi_g' * y / (phi_g' * phi_g + eps))
%
%     3. 重构残差（归一化）：
%           r_g = ||y - q_g* * phi_g||_2 / (||y||_2 + eps)
%
%     4. 综合距离：
%           D_g = lambda_shape * d_shape(g) + lambda_resid * r_g
%
%   输入：
%       F_obs_lin       - (M x 1) 线性功率观测
%       F_obs_shape_l1  - (M x 1) L1 归一化观测形状
%       SpatialFP_band  - 单频带指纹结构体，需含 F_lin (MxG), F_shape_l1 (MxG)
%       match_cfg       - [可选] 配置
%                           .lambda_shape (默认 0.7)
%                           .lambda_resid (默认 0.3)
%
%   输出：
%       D_vec       - (1 x G) 综合距离
%       d_shape_vec - (1 x G) 形状距离
%       resid_vec   - (1 x G) 重构残差
%       q_opt_vec   - (1 x G) 最优缩放系数

    % ---- 默认配置 ----
    if nargin < 4 || isempty(match_cfg)
        match_cfg = struct();
    end
    cfg = fill_shape_scale_defaults(match_cfg);

    F_lin   = SpatialFP_band.F_lin;         % M x G
    F_shape = SpatialFP_band.F_shape_l1;    % M x G

    % ---- 1) 形状距离 (L2) ----
    % d_shape(g) = ||s_obs - s_g||_2
    diff_shape  = F_obs_shape_l1 - F_shape;         % M x G
    d_shape_vec = sqrt(sum(diff_shape.^2, 1));       % 1 x G

    % ---- 2) 最优缩放系数 ----
    % q_g* = max(0, phi_g' * y / (phi_g' * phi_g + eps))
    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);        % 1 x G
    phi_dot_phi = sum(F_lin.^2, 1);                   % 1 x G
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + cfg.eps_val));  % 1 x G

    % ---- 3) 重构残差（归一化） ----
    % r_g = ||y - q_g* * phi_g||_2 / (||y||_2 + eps)
    residual    = F_obs_lin - q_opt_vec .* F_lin;     % M x G
    norm_y      = norm(F_obs_lin);
    resid_vec   = sqrt(sum(residual.^2, 1)) / (norm_y + cfg.eps_val);  % 1 x G

    % ---- 4) 综合距离 ----
    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;
end


%% ==================== 局部函数 ====================

function cfg = fill_shape_scale_defaults(cfg_in)
    cfg = cfg_in;

    if ~isfield(cfg, 'lambda_shape')
        cfg.lambda_shape = 0.7;
    end
    if ~isfield(cfg, 'lambda_resid')
        cfg.lambda_resid = 0.3;
    end
    if ~isfield(cfg, 'eps_val')
        cfg.eps_val = 1e-30;
    end
end
