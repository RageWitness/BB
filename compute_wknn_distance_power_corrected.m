function [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected( ...
    F_obs, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_POWER_CORRECTED  功率修正 WKNN 距离
%
%   [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected(
%       F_obs, SpatialFP_band)
%   [dist_vec, alpha_opt, F_obs_corr] = compute_wknn_distance_power_corrected(
%       F_obs, SpatialFP_band, match_cfg)
%
%   算法：
%     1. 对每个参考点 g，计算最优平移：
%           alpha_g* = mean(F_db(:,g) - F_obs)
%        然后裁剪到 [alpha_min, alpha_max]
%     2. 校正后的观测：F_obs_corr(:,g) = F_obs + alpha_g*
%     3. 距离计算（两种模式）：
%
%        shifted_euclidean（默认）：
%           d_g = ||F_obs + alpha_g*1 - F_db(:,g)||_2
%
%        shape_plus_shift：
%           d_g^2 = lambda_c * ||F_obs_c - F_db_c(:,g)||^2
%                 + lambda_mu * (alpha_g* - alpha_prior)^2
%
%        variance_weighted：
%           d_g^2 = sum_m (F_obs+alpha_g - F_db(m,g))^2 / (std_dBm(m,g)^2 + eps)
%
%     4. 若启用 Top-K AP：只用观测 RSS 最强的 K 个 AP 计算距离
%
%   输入：
%       F_obs          - (M x 1) RSS 观测 (dBm)
%       SpatialFP_band - 单频带指纹结构体
%       match_cfg      - [可选] 匹配配置，缺省自动填充默认值
%
%   输出：
%       dist_vec   - (1 x G) 与每个参考点的距离
%       alpha_opt  - (1 x G) 每个参考点的最优平移量 (dB)
%       F_obs_corr - (M x 1) 用全局中位 alpha 校正后的观测指纹

    % ---- 填充默认配置 ----
    if nargin < 3 || isempty(match_cfg)
        match_cfg = struct();
    end
    cfg = fill_match_defaults(match_cfg);

    M = numel(F_obs);
    G = size(SpatialFP_band.F_dBm, 2);

    % ---- Top-K AP 选取 ----
    if cfg.use_top_ap_only && cfg.top_ap_k < M
        [~, sort_idx] = sort(F_obs, 'descend');
        ap_mask = sort_idx(1:cfg.top_ap_k);
        F_obs_sub  = F_obs(ap_mask);
        F_db_sub   = SpatialFP_band.F_dBm(ap_mask, :);
        F_db_c_sub = SpatialFP_band.centered_dBm(ap_mask, :);
        if isfield(SpatialFP_band, 'std_dBm')
            Std_sub = SpatialFP_band.std_dBm(ap_mask, :);
        end
    else
        F_obs_sub  = F_obs;
        F_db_sub   = SpatialFP_band.F_dBm;
        F_db_c_sub = SpatialFP_band.centered_dBm;
        if isfield(SpatialFP_band, 'std_dBm')
            Std_sub = SpatialFP_band.std_dBm;
        end
    end

    % ---- 最优平移量 alpha_g* ----
    % alpha_g* = mean(F_db(:,g) - F_obs) 对选中的 AP
    mu_obs_sub = mean(F_obs_sub);
    mu_db_sub  = mean(F_db_sub, 1);          % 1 x G
    alpha_opt  = mu_db_sub - mu_obs_sub;     % 1 x G

    % 裁剪到约束区间
    alpha_opt = max(alpha_opt, cfg.alpha_min_dB);
    alpha_opt = min(alpha_opt, cfg.alpha_max_dB);

    % ---- 距离计算 ----
    switch cfg.mode
        case 'shifted_euclidean'
            % d_g = ||F_obs + alpha_g*1 - F_db(:,g)||
            % = ||(F_obs - F_db(:,g)) + alpha_g*1||
            diff_raw = F_db_sub - F_obs_sub;             % M_eff x G
            % diff_raw(:,g) = F_db(:,g) - F_obs
            % 校正后残差 = F_obs + alpha_g - F_db(:,g) = -(diff_raw - alpha_g)
            residual = diff_raw - alpha_opt;              % M_eff x G (广播)
            dist_vec = sqrt(sum(residual.^2, 1));         % 1 x G

        case 'shape_plus_shift'
            % 形状距离
            mu_obs_c = mean(F_obs_sub);
            F_obs_c  = F_obs_sub - mu_obs_c;             % M_eff x 1
            diff_shape = F_db_c_sub - F_obs_c;           % M_eff x G
            d_shape_sq = sum(diff_shape.^2, 1);           % 1 x G

            % 平移惩罚
            d_shift_sq = (alpha_opt - cfg.alpha_prior_dB).^2;  % 1 x G

            dist_vec = sqrt(cfg.lambda_c * d_shape_sq + ...
                            cfg.lambda_mu * d_shift_sq);

        case 'variance_weighted'
            % 方差加权平移欧氏距离
            % d_g^2 = sum_m (F_obs + alpha_g - F_db(m,g))^2 / (std(m,g)^2 + eps)
            if ~exist('Std_sub', 'var') || isempty(Std_sub)
                error('[M4] variance_weighted 模式需要 SpatialFP_band.std_dBm 字段');
            end
            diff_raw = F_db_sub - F_obs_sub;              % M_eff x G
            residual = diff_raw - alpha_opt;               % M_eff x G
            var_w = Std_sub.^2 + cfg.var_eps;              % M_eff x G
            dist_vec = sqrt(sum(residual.^2 ./ var_w, 1)); % 1 x G

        otherwise
            error('[M4] 未知 match_cfg.mode: %s', cfg.mode);
    end

    % ---- 输出校正后的观测指纹（用中位 alpha） ----
    alpha_median = median(alpha_opt);
    F_obs_corr = F_obs + alpha_median;
end


%% ==================== 局部函数 ====================

function cfg = fill_match_defaults(cfg_in)
% FILL_MATCH_DEFAULTS  填充匹配配置的默认值

    cfg = cfg_in;

    if ~isfield(cfg, 'mode')
        cfg.mode = 'shifted_euclidean';
    end

    % alpha 约束区间
    if ~isfield(cfg, 'alpha_min_dB')
        cfg.alpha_min_dB = -5;    % 默认：允许观测比库强 5 dB
    end
    if ~isfield(cfg, 'alpha_max_dB')
        cfg.alpha_max_dB = 25;    % 默认：允许观测比库弱 25 dB
    end

    % shape_plus_shift 模式参数
    if ~isfield(cfg, 'lambda_c')
        cfg.lambda_c = 1.0;
    end
    if ~isfield(cfg, 'lambda_mu')
        cfg.lambda_mu = 1.0;
    end
    if ~isfield(cfg, 'alpha_prior_dB')
        cfg.alpha_prior_dB = 0;
    end

    % variance_weighted 模式参数
    if ~isfield(cfg, 'var_eps')
        cfg.var_eps = 1.0;   % 防除零 + 正则化 (dB^2)
    end

    % Top-K AP
    if ~isfield(cfg, 'use_top_ap_only')
        cfg.use_top_ap_only = false;
    end
    if ~isfield(cfg, 'top_ap_k')
        cfg.top_ap_k = 8;
    end
end
