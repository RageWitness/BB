function [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected( ...
    F_obs_lin, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_POWER_CORRECTED  线性域功率修正 WKNN 距离
%
%   [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected(
%       F_obs_lin, SpatialFP_band)
%   [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected(
%       F_obs_lin, SpatialFP_band, match_cfg)
%
%   所有计算在线性功率域 (mW) 进行。
%
%   算法：
%     1. 对每个参考点 g，计算最优缩放因子（线性比例）：
%           alpha_g = mean(F_db_lin(:,g)) / mean(F_obs_lin)
%        裁剪到 [alpha_min, alpha_max]
%     2. 距离模式：
%
%        shifted_euclidean（默认）：
%           d_g = ||alpha_g * F_obs_lin - F_db_lin(:,g)||_2
%
%        shape_only：
%           归一化后比较形状：
%           d_g = ||F_obs_lin/mean(F_obs_lin) - F_db_lin(:,g)/mean(F_db_lin(:,g))||_2
%
%        variance_weighted：
%           d_g^2 = sum_m (alpha_g*F_obs(m) - F_db(m,g))^2 / (std_lin(m,g)^2 + eps)
%
%     3. 若启用 Top-K AP：只用观测 RSS 最强的 K 个 AP
%
%   输入：
%       F_obs_lin      - (M x 1) RSS 观测 (线性 mW)
%       SpatialFP_band - 单频带指纹结构体，需含 F_lin, mean_lin, centered_lin
%       match_cfg      - [可选] 匹配配置
%
%   输出：
%       dist_vec   - (1 x G) 距离
%       alpha_opt  - (1 x G) 每个参考点的最优缩放因子
%       F_obs_corr - (M x 1) 用中位 alpha 校正后的观测

    % ---- 默认配置 ----
    if nargin < 3 || isempty(match_cfg)
        match_cfg = struct();
    end
    cfg = fill_match_defaults(match_cfg);

    M = numel(F_obs_lin);

    % ---- Top-K AP 选取 ----
    if cfg.use_top_ap_only && cfg.top_ap_k < M
        [~, sort_idx] = sort(F_obs_lin, 'descend');
        ap_idx = sort_idx(1:cfg.top_ap_k);
        F_obs_sub = F_obs_lin(ap_idx);
        F_db_sub  = SpatialFP_band.F_lin(ap_idx, :);
        F_cn_sub  = SpatialFP_band.centered_lin(ap_idx, :);
        if isfield(SpatialFP_band, 'std_lin')
            Std_sub = SpatialFP_band.std_lin(ap_idx, :);
        end
    else
        F_obs_sub = F_obs_lin;
        F_db_sub  = SpatialFP_band.F_lin;
        F_cn_sub  = SpatialFP_band.centered_lin;
        if isfield(SpatialFP_band, 'std_lin')
            Std_sub = SpatialFP_band.std_lin;
        end
    end

    % ---- 最优缩放因子 alpha_g ----
    % alpha_g = mean(F_db(:,g)) / mean(F_obs)
    mu_obs = mean(F_obs_sub);
    mu_obs = max(mu_obs, 1e-30);  % 防除零
    mu_db  = mean(F_db_sub, 1);   % 1 x G
    alpha_opt = mu_db / mu_obs;   % 1 x G（线性比例）

    % 裁剪
    alpha_opt = max(alpha_opt, cfg.alpha_min);
    alpha_opt = min(alpha_opt, cfg.alpha_max);

    % ---- 距离计算 ----
    switch cfg.mode
        case 'shifted_euclidean'
            % d_g = ||alpha_g * F_obs - F_db(:,g)||
            residual = alpha_opt .* F_obs_sub - F_db_sub;  % M_eff x G
            dist_vec = sqrt(sum(residual.^2, 1));           % 1 x G

        case 'shape_only'
            % 归一化形状比较
            F_obs_n = F_obs_sub / mu_obs;                   % M_eff x 1
            diff_shape = F_cn_sub - F_obs_n;                % M_eff x G
            dist_vec = sqrt(sum(diff_shape.^2, 1));         % 1 x G

        case 'variance_weighted'
            if ~exist('Std_sub', 'var') || isempty(Std_sub)
                error('[M4] variance_weighted 需要 SpatialFP_band.std_lin');
            end
            residual = alpha_opt .* F_obs_sub - F_db_sub;
            var_w = Std_sub.^2 + cfg.var_eps;
            dist_vec = sqrt(sum(residual.^2 ./ var_w, 1));

        otherwise
            error('[M4] 未知 match_cfg.mode: %s', cfg.mode);
    end

    % ---- 校正后观测 ----
    alpha_median = median(alpha_opt);
    F_obs_corr = alpha_median * F_obs_lin;
end


%% ==================== 局部函数 ====================

function cfg = fill_match_defaults(cfg_in)
    cfg = cfg_in;

    if ~isfield(cfg, 'mode')
        cfg.mode = 'shifted_euclidean';
    end

    % alpha 约束区间（线性比例）
    % 库 ref_power = 100 dBm, 目标可能 50-140 dBm
    %   alpha = P_ref_lin / P_tx_lin = 10^((P_ref - P_tx)/10)
    %   P_tx = 140 → alpha = 10^(-40/10) = 0.0001
    %   P_tx = 50  → alpha = 10^(50/10)  = 100000
    %   宽松默认：[1e-5, 1e6]
    if ~isfield(cfg, 'alpha_min')
        cfg.alpha_min = 1e-5;
    end
    if ~isfield(cfg, 'alpha_max')
        cfg.alpha_max = 1e6;
    end

    % variance_weighted 防除零
    if ~isfield(cfg, 'var_eps')
        cfg.var_eps = 1e-20;
    end

    % Top-K AP
    if ~isfield(cfg, 'use_top_ap_only')
        cfg.use_top_ap_only = false;
    end
    if ~isfield(cfg, 'top_ap_k')
        cfg.top_ap_k = 8;
    end
end
