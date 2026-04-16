function [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = compute_wknn_distance_shape_scale( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_SHAPE_SCALE  基于 shape + scale 的 WKNN 距离
%
%   支持模式（match_cfg.distance_mode）：
%     'shape_scale'         - 原始模式
%     'shape_scale_prob'    - 概率 shape + 全局 AP 权重
%     'shape_scale_hn'      - hard-negative 概率增强模式
%     'shape_scale_masked'  - 可靠 AP 掩码模式

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

        case 'shape_scale_masked'
            [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = ...
                compute_masked_mode(F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg);

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

    residual  = F_obs_lin - q_opt_vec .* F_lin;
    norm_y    = norm(F_obs_lin);
    resid_vec = sqrt(sum(residual.^2, 1)) / (norm_y + cfg.eps_val);

    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;
end


function [D_vec, d_shape_vec, resid_vec, q_opt_vec, mask_info] = compute_masked_mode( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg)
% shape_scale_masked：
%   1) 可靠 AP 掩码 v（由 excess over noise floor 计算）
%   2) 掩码 shape 距离
%   3) q* 只在掩码 AP 计算
%   4) residual 只在掩码 AP 计算

    F_lin = SpatialFP_band.F_lin;
    [M, ~] = size(F_lin);

    if numel(F_obs_lin) ~= M
        error('[M4] shape_scale_masked: F_obs_lin size does not match F_lin');
    end

    mask_info = struct();
    mask_info.masked_used = false;
    mask_info.fallback_used = false;
    mask_info.fallback_reason = '';

    % 关闭 masked 时，直接走旧逻辑
    if ~cfg.masked_enable
        [D_vec, d_shape_vec, resid_vec, q_opt_vec] = ...
            compute_original_mode(F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg);
        mask_info.fallback_used = true;
        mask_info.fallback_reason = 'masked_disabled';
        mask_info.valid_ap_mask = true(M, 1);
        mask_info.valid_ap_weight = ones(M, 1);
        mask_info.n_valid_ap = M;
        mask_info.d_shape_mask_vec = d_shape_vec;
        mask_info.resid_mask_vec = resid_vec;
        return;
    end

    % 观测 dBm
    if isfield(cfg, 'F_obs_dBm') && ~isempty(cfg.F_obs_dBm)
        F_obs_dBm = cfg.F_obs_dBm(:);
    else
        F_obs_dBm = 10 * log10(max(F_obs_lin(:), cfg.masked_eps));
    end

    % 无噪声底，回退旧逻辑（并记录原因）
    if ~isfield(cfg, 'noise_floor_dBm') || isempty(cfg.noise_floor_dBm)
        [D_vec, d_shape_vec, resid_vec, q_opt_vec] = ...
            compute_original_mode(F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg);
        mask_info.fallback_used = true;
        mask_info.fallback_reason = 'missing_noise_floor';
        mask_info.valid_ap_mask = true(M, 1);
        mask_info.valid_ap_weight = ones(M, 1);
        mask_info.n_valid_ap = M;
        mask_info.d_shape_mask_vec = d_shape_vec;
        mask_info.resid_mask_vec = resid_vec;
        return;
    end

    noise_floor = cfg.noise_floor_dBm(:);
    if numel(noise_floor) == 1
        noise_floor = repmat(noise_floor, M, 1);
    elseif numel(noise_floor) ~= M
        error('[M4] shape_scale_masked: noise_floor_dBm length must be 1 or M');
    end

    if cfg.tau_high_dB <= cfg.tau_low_dB
        cfg.tau_high_dB = cfg.tau_low_dB + 1e-6;
    end

    % reliable mask
    excess_dB = F_obs_dBm - noise_floor;
    v = zeros(M, 1);
    idx_hi  = excess_dB >= cfg.tau_high_dB;
    idx_mid = excess_dB > cfg.tau_low_dB & excess_dB < cfg.tau_high_dB;
    v(idx_hi) = 1;
    v(idx_mid) = (excess_dB(idx_mid) - cfg.tau_low_dB) / ...
                 (cfg.tau_high_dB - cfg.tau_low_dB + cfg.masked_eps);
    v = max(0, min(1, v));

    valid_mask = excess_dB >= cfg.tau_low_dB;
    n_valid_ap = sum(valid_mask);

    % 全部无效时回退
    if sum(v) <= cfg.masked_eps
        [D_vec, d_shape_vec, resid_vec, q_opt_vec] = ...
            compute_original_mode(F_obs_lin, F_obs_shape_l1, SpatialFP_band, cfg);
        mask_info.fallback_used = true;
        mask_info.fallback_reason = 'no_valid_ap';
        mask_info.valid_ap_mask = valid_mask;
        mask_info.valid_ap_weight = v;
        mask_info.n_valid_ap = n_valid_ap;
        mask_info.excess_dB = excess_dB;
        mask_info.noise_floor_dBm = noise_floor;
        mask_info.d_shape_mask_vec = d_shape_vec;
        mask_info.resid_mask_vec = resid_vec;
        return;
    end

    % --- d_shape_mask ---
    y_mask = v .* F_obs_lin;                           % M x 1
    s_obs_mask = y_mask / (sum(y_mask) + cfg.masked_eps);
    phi_mask = F_lin .* v;                             % M x G
    s_g_mask = phi_mask ./ (sum(phi_mask, 1) + cfg.masked_eps);
    d_shape_vec = sqrt(sum((s_obs_mask - s_g_mask).^2, 1));

    % --- q*_mask ---
    phi_dot_y   = sum(phi_mask .* y_mask, 1);
    phi_dot_phi = sum(phi_mask.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + cfg.masked_eps));

    % --- residual_mask ---
    residual  = v .* (F_obs_lin - q_opt_vec .* F_lin);
    norm_y_w  = norm(y_mask);
    resid_vec = sqrt(sum(residual.^2, 1)) / (norm_y_w + cfg.masked_eps);

    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;

    mask_info.masked_used = true;
    mask_info.valid_ap_mask = valid_mask;
    mask_info.valid_ap_weight = v;
    mask_info.n_valid_ap = n_valid_ap;
    mask_info.excess_dB = excess_dB;
    mask_info.noise_floor_dBm = noise_floor;
    mask_info.d_shape_mask_vec = d_shape_vec;
    mask_info.resid_mask_vec = resid_vec;
end


function [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = compute_prob_mode( ...
    F_obs_lin, SpatialFP_band, cfg)
% 概率增强 shape_scale_prob 模式

    F_lin     = SpatialFP_band.F_lin;
    mu_shape  = SpatialFP_band.F_shape_prob_mu;
    var_shape = SpatialFP_band.F_shape_prob_var;
    w_ap      = SpatialFP_band.ap_weight_global;
    std_lin   = SpatialFP_band.std_lin;

    M = size(F_lin, 1);
    eps_p = cfg.prob_eps;

    y_w = w_ap .* F_obs_lin;
    s_obs_w = y_w / (sum(y_w) + eps_p);

    diff_mu = s_obs_w - mu_shape;
    d_shape_vec = sqrt(sum(diff_mu.^2 ./ (var_shape + eps_p), 1));

    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);
    phi_dot_phi = sum(F_lin.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + eps_p));

    residual  = F_obs_lin - q_opt_vec .* F_lin;
    resid_vec = sqrt(sum(residual.^2 ./ (std_lin.^2 + eps_p), 1) / M);

    D_vec = cfg.lambda_shape * d_shape_vec + cfg.lambda_resid * resid_vec;

    prob_info.s_obs_weighted   = s_obs_w;
    prob_info.d_shape_prob_vec = d_shape_vec;
    prob_info.resid_prob_vec   = resid_vec;
    prob_info.ap_weight_used   = w_ap;
end


function [D_vec, d_shape_vec, resid_vec, q_opt_vec, hn_info] = compute_hn_mode( ...
    F_obs_lin, SpatialFP_band, cfg)
% 局部 hard-negative 增强模式 shape_scale_hn

    F_lin     = SpatialFP_band.F_lin;              % M x G_sub
    mu_shape  = SpatialFP_band.F_shape_prob_mu;    % M x G_sub
    var_shape = SpatialFP_band.F_shape_prob_var;   % M x G_sub
    w_hn      = SpatialFP_band.ap_weight_hn;       % M x G_sub
    std_lin   = SpatialFP_band.std_lin;            % M x G_sub
    hn_idx    = SpatialFP_band.hardneg_idx;        % G_sub x top_conf

    % 全库 mu/var（hardneg_idx 指向全库索引）
    if isfield(SpatialFP_band, 'full_F_shape_prob_mu')
        mu_full  = SpatialFP_band.full_F_shape_prob_mu;
        var_full = SpatialFP_band.full_F_shape_prob_var;
    else
        mu_full  = mu_shape;
        var_full = var_shape;
    end

    [M, G_sub] = size(F_lin);
    eps_p = cfg.prob_eps;
    delta = cfg.hn_delta_margin;
    lam_s = cfg.lambda_shape;
    lam_r = cfg.lambda_resid;
    lam_h = cfg.lambda_hn;

    d_shape_vec   = zeros(1, G_sub);
    d_neg_min_vec = inf(1, G_sub);
    p_hn_vec      = zeros(1, G_sub);

    for g = 1:G_sub
        w_g = w_hn(:, g);
        y_w = w_g .* F_obs_lin;
        s_obs_hn = y_w / (sum(y_w) + eps_p);

        diff_g = s_obs_hn - mu_shape(:, g);
        d_shape_hn_g = sqrt(sum(diff_g.^2 ./ (var_shape(:, g) + eps_p)));
        d_shape_vec(g) = d_shape_hn_g;

        hn_g = hn_idx(g, :);
        hn_g = hn_g(hn_g > 0);  % 去零填充
        if ~isempty(hn_g)
            mu_H  = mu_full(:, hn_g);
            var_H = var_full(:, hn_g);
            diff_H = s_obs_hn - mu_H;
            d_neg_h = sqrt(sum(diff_H.^2 ./ (var_H + eps_p), 1));
            d_neg_min_vec(g) = min(d_neg_h);
            p_hn_vec(g) = max(0, delta - (d_neg_min_vec(g) - d_shape_hn_g));
        end
    end

    phi_dot_y   = sum(F_lin .* F_obs_lin, 1);
    phi_dot_phi = sum(F_lin.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + eps_p));

    residual  = F_obs_lin - q_opt_vec .* F_lin;
    resid_vec = sqrt(sum(residual.^2 ./ (std_lin.^2 + eps_p), 1) / M);

    D_vec = lam_s * d_shape_vec + lam_r * resid_vec + lam_h * p_hn_vec;

    hn_info.d_shape_hn_vec = d_shape_vec;
    hn_info.resid_prob_vec = resid_vec;
    hn_info.d_neg_min_vec  = d_neg_min_vec;
    hn_info.p_hn_vec       = p_hn_vec;
end


function cfg = fill_shape_scale_defaults(cfg_in)
% FILL_SHAPE_SCALE_DEFAULTS  填充匹配参数默认值
    cfg = cfg_in;

    if ~isfield(cfg, 'distance_mode')
        cfg.distance_mode = 'shape_scale';
    end

    if ~isfield(cfg, 'lambda_shape')
        switch cfg.distance_mode
            case 'shape_scale_hn'
                cfg.lambda_shape = 0.50;
            case 'shape_scale_prob'
                cfg.lambda_shape = 0.60;
            otherwise
                cfg.lambda_shape = 0.70;
        end
    end

    if ~isfield(cfg, 'lambda_resid')
        switch cfg.distance_mode
            case 'shape_scale_hn'
                cfg.lambda_resid = 0.25;
            case 'shape_scale_prob'
                cfg.lambda_resid = 0.40;
            otherwise
                cfg.lambda_resid = 0.30;
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

    % masked 参数：支持 nested + flat 两种传参
    cfg.masked_enable = get_field_with_fallback(cfg, {'masked_shape', 'enable'}, true);
    cfg.tau_low_dB    = get_field_with_fallback(cfg, {'masked_shape', 'tau_low_dB'}, 3);
    cfg.tau_high_dB   = get_field_with_fallback(cfg, {'masked_shape', 'tau_high_dB'}, 10);
    cfg.masked_eps    = get_field_with_fallback(cfg, {'masked_shape', 'eps'}, 1e-12);

    if isfield(cfg_in, 'masked_enable'); cfg.masked_enable = cfg_in.masked_enable; end
    if isfield(cfg_in, 'tau_low_dB');    cfg.tau_low_dB = cfg_in.tau_low_dB; end
    if isfield(cfg_in, 'tau_high_dB');   cfg.tau_high_dB = cfg_in.tau_high_dB; end
    if isfield(cfg_in, 'masked_eps');    cfg.masked_eps = cfg_in.masked_eps; end
end


function v = get_field_with_fallback(s, path, def)
% GET_FIELD_WITH_FALLBACK  安全读取嵌套字段
    cur = s;
    for i = 1:numel(path)
        fn = path{i};
        if isstruct(cur) && isfield(cur, fn)
            cur = cur.(fn);
        else
            v = def;
            return;
        end
    end
    v = cur;
end
