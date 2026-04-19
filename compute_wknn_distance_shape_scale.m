function [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = compute_wknn_distance_shape_scale( ...
    F_obs_lin, F_obs_shape_l1, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_SHAPE_SCALE  基于 shape + scale 的 WKNN 距离
%
%   支持模式（match_cfg.distance_mode）：
%     'shape_scale'         - 原始 shape + residual
%     'shape_scale_masked'  - 可靠 AP 掩码下的 shape + residual

    if nargin < 4 || isempty(match_cfg)
        match_cfg = struct();
    end
    cfg = fill_shape_scale_defaults(match_cfg);
    prob_info = struct();

    switch cfg.distance_mode
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

    if isfield(cfg, 'F_obs_dBm') && ~isempty(cfg.F_obs_dBm)
        F_obs_dBm = cfg.F_obs_dBm(:);
    else
        F_obs_dBm = 10 * log10(max(F_obs_lin(:), cfg.masked_eps));
    end

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

    y_mask = v .* F_obs_lin;
    s_obs_mask = y_mask / (sum(y_mask) + cfg.masked_eps);
    phi_mask = F_lin .* v;
    s_g_mask = phi_mask ./ (sum(phi_mask, 1) + cfg.masked_eps);
    d_shape_vec = sqrt(sum((s_obs_mask - s_g_mask).^2, 1));

    phi_dot_y   = sum(phi_mask .* y_mask, 1);
    phi_dot_phi = sum(phi_mask.^2, 1);
    q_opt_vec   = max(0, phi_dot_y ./ (phi_dot_phi + cfg.masked_eps));

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


function cfg = fill_shape_scale_defaults(cfg_in)
    cfg = cfg_in;
    if ~isfield(cfg, 'distance_mode'), cfg.distance_mode = 'shape_scale'; end
    if ~isfield(cfg, 'lambda_shape'),  cfg.lambda_shape  = 0.70; end
    if ~isfield(cfg, 'lambda_resid'),  cfg.lambda_resid  = 0.30; end
    if ~isfield(cfg, 'eps_val'),       cfg.eps_val       = 1e-30; end

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
    cur = s;
    for i = 1:numel(path)
        fn = path{i};
        if isstruct(cur) && isfield(cur, fn)
            cur = cur.(fn);
        else
            v = def; return;
        end
    end
    v = cur;
end
