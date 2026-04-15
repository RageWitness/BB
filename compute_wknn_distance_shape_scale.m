function [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = compute_wknn_distance_shape_scale( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_SHAPE_SCALE  基于形状+缩放的 WKNN 距离
%
%   支持三种距离模式（通过 match_cfg.distance_mode 切换）：
%
%   'shape_scale'（默认，原始模式）
%   'shape_scale_prob'（全局 AP 权重 + 概率 shape）
%   'shape_scale_hn'（局部 hard-negative AP 权重 + 排斥项）
%
%   输入：
%       F_obs_lin       - (M x 1) 线性功率观测
%       F_obs_shape_l1  - (M x 1) L1 归一化观测形状
%       SpatialFP_band  - 单频带指纹结构体
%       match_cfg       - [可选] 配置
%
%   输出：
%       D_vec       - (1 x G) 综合距离
%       d_shape_vec - (1 x G) 形状距离
%       resid_vec   - (1 x G) 重构残差
%       q_opt_vec   - (1 x G) 最优缩放系数
%       prob_info   - [可选] 额外诊断

    if nargin < 4 || isempty(match_cfg)
        match_cfg = struct();
    end
    cfg = fill_shape_scale_defaults(match_cfg);

    prob_info = struct();

    switch cfg.distance_mode
        case 'shape_scale_hn'
            [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = ...
                compute_hn_mode(F_obs_lin, SpatialFP_band, cfg);
        case 'shape_scale_prob'
            [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = ...
                compute_prob_mode(F_obs_lin, SpatialFP_band, cfg);
        otherwise
            [D_vec, d_shape_vec, resid_vec, q_opt_vec] = ...
                compute_original_mode(F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg);
    end
end


%% ==================== 局部函数 ====================

function [D_vec, d_shape_vec, resid_vec, q_opt_vec] = compute_original_mode( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg)
% 原始 shape_scale 模式

    F_lin   = SpatialFP_band.F_lin;
    F_shape = SpatialFP_band.F_shape_l1;

    diff_shape  = F_obs_shape_l1 - F_shape;
    d_shape_vec = sqrt(sum(diff_shape.^2, 1));

    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);
    phi_dot_phi = sum(F_lin.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + cfg.eps_val));

    residual    = F_obs_lin - q_opt_vec .* F_lin;
    norm_y      = norm(F_obs_lin);
    resid_vec   = sqrt(sum(residual.^2, 1)) / (norm_y + cfg.eps_val);

    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;
end


function [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = compute_prob_mode( ...
    F_obs_lin, SpatialFP_band, cfg)
% 概率增强 shape_scale_prob 模式

    F_lin    = SpatialFP_band.F_lin;
    mu_shape = SpatialFP_band.F_shape_prob_mu;
    var_shape = SpatialFP_band.F_shape_prob_var;
    w_ap     = SpatialFP_band.ap_weight_global;
    std_lin  = SpatialFP_band.std_lin;

    M = size(F_lin, 1);
    eps_p = cfg.prob_eps;

    y_w = w_ap .* F_obs_lin;
    s_obs_w = y_w / (sum(y_w) + eps_p);

    diff_mu     = s_obs_w - mu_shape;
    d_shape_vec = sqrt(sum(diff_mu.^2 ./ (var_shape + eps_p), 1));

    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);
    phi_dot_phi = sum(F_lin.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + eps_p));

    residual    = F_obs_lin - q_opt_vec .* F_lin;
    resid_vec   = sqrt(sum(residual.^2 ./ (std_lin.^2 + eps_p), 1) / M);

    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;

    prob_info.s_obs_weighted   = s_obs_w;
    prob_info.d_shape_prob_vec = d_shape_vec;
    prob_info.resid_prob_vec   = resid_vec;
    prob_info.ap_weight_used   = w_ap;
end


function [D_vec, d_shape_vec, resid_vec, q_opt_vec, hn_info] = compute_hn_mode( ...
    F_obs_lin, SpatialFP_band, cfg)
% 局部 hard-negative 增强模式 shape_scale_hn
%
%   对每个候选 g：
%     1. 读局部权重 w_hn(:,g) → 构造 s_obs_hn^(g)
%     2. d_shape_hn(g) = sqrt(sum (s_obs_hn - mu)^2 / (var + eps))
%     3. q_g^* 保持原定义
%     4. r_prob(g) = sqrt((1/M) sum (y - q*phi)^2 / (std^2 + eps))
%     5. d_neg_min(g): 观测对 g 的 hard-neg 的最小概率 shape 距离
%     6. p_hn(g) = max(0, delta - (d_neg_min - d_shape_hn))
%     7. D_hn(g) = lam_shape*d_shape_hn + lam_resid*r_prob + lam_hn*p_hn

    F_lin    = SpatialFP_band.F_lin;              % M x G_sub
    mu_shape = SpatialFP_band.F_shape_prob_mu;    % M x G_sub
    var_shape = SpatialFP_band.F_shape_prob_var;  % M x G_sub
    w_hn     = SpatialFP_band.ap_weight_hn;       % M x G_sub
    std_lin  = SpatialFP_band.std_lin;            % M x G_sub
    hn_idx   = SpatialFP_band.hardneg_idx;        % G_sub x top_conf

    % 全库 mu/var（hardneg_idx 指向全库索引，需用全库查找）
    if isfield(SpatialFP_band, 'full_F_shape_prob_mu')
        mu_full  = SpatialFP_band.full_F_shape_prob_mu;   % M x G_full
        var_full = SpatialFP_band.full_F_shape_prob_var;   % M x G_full
    else
        mu_full  = mu_shape;   % 非子集调用时，就是全库自身
        var_full = var_shape;
    end

    [M, G_sub] = size(F_lin);
    eps_p    = cfg.prob_eps;
    delta    = cfg.hn_delta_margin;
    lam_s    = cfg.lambda_shape;
    lam_r    = cfg.lambda_resid;
    lam_h    = cfg.lambda_hn;

    % 预分配
    d_shape_vec = zeros(1, G_sub);
    d_neg_min_vec = inf(1, G_sub);
    p_hn_vec    = zeros(1, G_sub);

    for g = 1:G_sub
        w_g = w_hn(:, g);                         % M x 1  局部权重

        % ---- 候选相关加权观测 shape ----
        y_w = w_g .* F_obs_lin;                    % M x 1
        s_obs_hn = y_w / (sum(y_w) + eps_p);      % M x 1

        % ---- d_shape_hn(g) ----
        diff_g = s_obs_hn - mu_shape(:, g);        % M x 1
        d_shape_hn_g = sqrt(sum(diff_g.^2 ./ (var_shape(:, g) + eps_p)));
        d_shape_vec(g) = d_shape_hn_g;

        % ---- hard-negative 排斥项 ----
        hn_g = hn_idx(g, :);                       % 1 x top_conf
        hn_g = hn_g(hn_g > 0);                     % 过滤零填充

        if ~isempty(hn_g)
            % hardneg_idx 指向全库索引，用全库 mu/var 查找
            mu_H = mu_full(:, hn_g);              % M x n_hn
            var_H = var_full(:, hn_g);            % M x n_hn

            diff_H = s_obs_hn - mu_H;             % M x n_hn
            d_neg_h = sqrt(sum(diff_H.^2 ./ (var_H + eps_p), 1));  % 1 x n_hn
            d_neg_min_vec(g) = min(d_neg_h);

            p_hn_vec(g) = max(0, delta - (d_neg_min_vec(g) - d_shape_hn_g));
        end
        % 若 hn_g 为空，d_neg_min = inf, p_hn = 0（已初始化）
    end

    % ---- q_g^* 保持原定义（向量化） ----
    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);
    phi_dot_phi = sum(F_lin.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + eps_p));

    % ---- 概率残差 r_prob (向量化) ----
    residual  = F_obs_lin - q_opt_vec .* F_lin;
    resid_vec = sqrt(sum(residual.^2 ./ (std_lin.^2 + eps_p), 1) / M);

    % ---- D_hn 综合距离 ----
    D_vec = lam_s * d_shape_vec + lam_r * resid_vec + lam_h * p_hn_vec;

    % ---- 诊断输出 ----
    hn_info.d_shape_hn_vec = d_shape_vec;
    hn_info.resid_prob_vec = resid_vec;
    hn_info.d_neg_min_vec  = d_neg_min_vec;
    hn_info.p_hn_vec       = p_hn_vec;
end


function cfg = fill_shape_scale_defaults(cfg_in)
    cfg = cfg_in;

    if ~isfield(cfg, 'distance_mode')
        cfg.distance_mode = 'shape_scale';
    end
    if ~isfield(cfg, 'lambda_shape')
        switch cfg.distance_mode
            case 'shape_scale_hn'
                cfg.lambda_shape = 0.50;
            case 'shape_scale_prob'
                cfg.lambda_shape = 0.6;
            otherwise
                cfg.lambda_shape = 0.7;
        end
    end
    if ~isfield(cfg, 'lambda_resid')
        switch cfg.distance_mode
            case 'shape_scale_hn'
                cfg.lambda_resid = 0.25;
            case 'shape_scale_prob'
                cfg.lambda_resid = 0.4;
            otherwise
                cfg.lambda_resid = 0.3;
        end
    end
    if ~isfield(cfg, 'lambda_hn')
        cfg.lambda_hn = 0.25;
    end
    if ~isfield(cfg, 'hn_delta_margin')
        cfg.hn_delta_margin = 0.15;
    end
    if ~isfield(cfg, 'eps_val')
        cfg.eps_val = 1e-30;
    end
    if ~isfield(cfg, 'prob_eps')
        cfg.prob_eps = 1e-12;
    end
end
