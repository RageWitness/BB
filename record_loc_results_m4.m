function loc_result = record_loc_results_m4(event, est_pos_xy, F_obs, ...
    neighbor_idx, neighbor_dist, weights, FrameStates, Config, extra_info)
% RECORD_LOC_RESULTS_M4  记录单个事件的定位结果

    loc_result = struct();

    % 基本信息
    loc_result.event_id     = event.event_id;
    loc_result.band_id      = event.band_id;
    loc_result.type_hat     = event.type_hat;
    loc_result.route_action = event.route_action;
    loc_result.time_range   = event.time_range;

    if isfield(event, 'label')
        loc_result.label = event.label;
    else
        loc_result.label = [];
    end
    if isfield(event, 'source_uid')
        loc_result.source_uid = event.source_uid;
    else
        loc_result.source_uid = '';
    end
    if isfield(event, 'linked_template_key')
        loc_result.linked_template_key = event.linked_template_key;
    else
        loc_result.linked_template_key = '';
    end

    % --- prior / metadata 透传 ---
    if isfield(event, 'location_prior') && isstruct(event.location_prior)
        loc_result.location_prior = event.location_prior;
    else
        loc_result.location_prior = struct('type', 'none', 'value', []);
    end
    if isfield(event, 'power_prior') && isstruct(event.power_prior)
        loc_result.power_prior = event.power_prior;
    else
        loc_result.power_prior = struct('type', 'none', 'value', []);
    end
    if isfield(event, 'metadata')
        loc_result.metadata = event.metadata;
    else
        loc_result.metadata = struct();
    end

    % 便利字段
    if isfield(event, 'source_type_name')
        loc_result.source_type_name = event.source_type_name;
    else
        loc_result.source_type_name = event.type_hat;
    end
    loc_result.location_prior_type = loc_result.location_prior.type;
    loc_result.power_prior_type    = loc_result.power_prior.type;
    if isfield(event, 'instance_id')
        loc_result.instance_id = event.instance_id;
    else
        loc_result.instance_id = 0;
    end
    if isfield(event, 'is_calibration_source')
        loc_result.is_calibration_source = event.is_calibration_source;
    else
        loc_result.is_calibration_source = false;
    end
    if isfield(event, 'is_target_source')
        loc_result.is_target_source = event.is_target_source;
    else
        loc_result.is_target_source = false;
    end
    if isfield(event, 'is_opportunistic_source')
        loc_result.is_opportunistic_source = event.is_opportunistic_source;
    else
        loc_result.is_opportunistic_source = false;
    end

    % 观测与估计
    loc_result.obs_fp_lin   = F_obs;
    loc_result.obs_fp_dBm   = 10 * log10(max(F_obs, 1e-30));
    loc_result.est_pos_xy   = est_pos_xy;

    % 近邻信息
    loc_result.neighbor_idx  = neighbor_idx;
    loc_result.neighbor_dist = neighbor_dist;
    loc_result.weights       = weights;
    loc_result.match_score   = 1 / (min(neighbor_dist) + 1e-10);

    % 真值与误差
    true_pos = extract_true_position(event, FrameStates);
    loc_result.true_pos_xy = true_pos;
    if ~isempty(true_pos) && ~all(true_pos == 0)
        loc_result.loc_error = norm(est_pos_xy - true_pos, 2);
    else
        loc_result.loc_error = NaN;
    end

    % 噪声底
    if isfield(Config, 'm1') && isfield(Config.m1, 'noise') && ...
            isfield(Config.m1.noise, 'n0_dBmHz')
        loc_result.n0 = Config.m1.noise.n0_dBmHz;
    else
        loc_result.n0 = -174;
    end

    % 默认诊断输出
    loc_result.q_opt_neighbors   = [];
    loc_result.d_shape_neighbors = [];
    loc_result.resid_neighbors   = [];
    loc_result.distance_mode     = 'shape_scale';
    loc_result.preselect_idx     = [];

    % 新增：masked 诊断
    loc_result.valid_ap_mask    = [];
    loc_result.valid_ap_weight  = [];
    loc_result.n_valid_ap       = [];
    loc_result.d_shape_mask_vec = [];
    loc_result.resid_mask_vec   = [];

    % prob/hn 诊断
    loc_result.prob_s_obs_weighted    = [];
    loc_result.prob_ap_weight_used    = [];
    loc_result.prob_d_shape_neighbors = [];
    loc_result.prob_resid_neighbors   = [];
    loc_result.hn_d_shape_neighbors   = [];
    loc_result.hn_resid_prob_neighbors = [];
    loc_result.hn_d_neg_min_neighbors = [];
    loc_result.hn_p_hn_neighbors      = [];

    % area screening 诊断
    loc_result.area_screening_enabled         = false;
    loc_result.area_screening_applied         = false;
    loc_result.area_screening_fallback_used   = false;
    loc_result.area_screening_fallback_reason = '';
    loc_result.area_screening_center_xy       = [NaN, NaN];
    loc_result.area_screening_radius_m        = NaN;
    loc_result.area_screening_pair_dist12_m   = NaN;
    loc_result.area_screening_candidate_xy    = [];
    loc_result.area_screening_candidate_dist_to_center_m = [];
    loc_result.area_screening_pre_idx_shape   = [];
    loc_result.area_screening_pre_dist_shape  = [];
    loc_result.area_screening_valid_idx       = [];
    loc_result.area_screening_rejected_idx    = [];
    loc_result.area_screening_valid_count     = 0;

    if nargin < 9 || isempty(extra_info)
        return;
    end

    % 基础匹配诊断
    if isfield(extra_info, 'q_opt_vec') && ~isempty(extra_info.q_opt_vec)
        loc_result.q_opt_neighbors = extra_info.q_opt_vec(neighbor_idx);
    end
    if isfield(extra_info, 'd_shape_vec') && ~isempty(extra_info.d_shape_vec)
        loc_result.d_shape_neighbors = extra_info.d_shape_vec(neighbor_idx);
    end
    if isfield(extra_info, 'resid_vec') && ~isempty(extra_info.resid_vec)
        loc_result.resid_neighbors = extra_info.resid_vec(neighbor_idx);
    end
    if isfield(extra_info, 'distance_mode')
        loc_result.distance_mode = extra_info.distance_mode;
    end
    if isfield(extra_info, 'preselect_idx')
        loc_result.preselect_idx = extra_info.preselect_idx;
    end
    if isfield(extra_info, 'valid_ap_mask')
        loc_result.valid_ap_mask = logical(extra_info.valid_ap_mask(:)');
    end
    if isfield(extra_info, 'n_valid_ap')
        loc_result.n_valid_ap = extra_info.n_valid_ap;
    end

    % prob/hn/masked 细项诊断
    if isfield(extra_info, 'prob_info') && ~isempty(fieldnames(extra_info.prob_info))
        pi = extra_info.prob_info;

        % prob
        if isfield(pi, 's_obs_weighted')
            loc_result.prob_s_obs_weighted = pi.s_obs_weighted;
        end
        if isfield(pi, 'ap_weight_used')
            loc_result.prob_ap_weight_used = pi.ap_weight_used;
        end
        if isfield(pi, 'd_shape_prob_vec')
            loc_result.prob_d_shape_neighbors = pi.d_shape_prob_vec(neighbor_idx);
        end
        if isfield(pi, 'resid_prob_vec')
            loc_result.prob_resid_neighbors = pi.resid_prob_vec(neighbor_idx);
        end

        % hn
        if isfield(pi, 'd_shape_hn_vec')
            loc_result.hn_d_shape_neighbors = pi.d_shape_hn_vec(neighbor_idx);
        end
        if isfield(pi, 'resid_prob_vec')
            loc_result.hn_resid_prob_neighbors = pi.resid_prob_vec(neighbor_idx);
        end
        if isfield(pi, 'd_neg_min_vec')
            loc_result.hn_d_neg_min_neighbors = pi.d_neg_min_vec(neighbor_idx);
        end
        if isfield(pi, 'p_hn_vec')
            loc_result.hn_p_hn_neighbors = pi.p_hn_vec(neighbor_idx);
        end

        % masked
        if isfield(pi, 'valid_ap_mask')
            loc_result.valid_ap_mask = logical(pi.valid_ap_mask(:)');
        end
        if isfield(pi, 'valid_ap_weight')
            loc_result.valid_ap_weight = pi.valid_ap_weight(:)';
        end
        if isfield(pi, 'n_valid_ap')
            loc_result.n_valid_ap = pi.n_valid_ap;
        elseif ~isempty(loc_result.valid_ap_mask)
            loc_result.n_valid_ap = sum(loc_result.valid_ap_mask);
        end
        if isfield(pi, 'd_shape_mask_vec')
            loc_result.d_shape_mask_vec = pi.d_shape_mask_vec;
        end
        if isfield(pi, 'resid_mask_vec')
            loc_result.resid_mask_vec = pi.resid_mask_vec;
        end
    end

    if isempty(loc_result.valid_ap_weight) && ~isempty(loc_result.valid_ap_mask)
        loc_result.valid_ap_weight = double(loc_result.valid_ap_mask);
    end

    % area screening 诊断（所有坐标/距离均为物理空间量）
    if isfield(extra_info, 'area_scr')
        as = extra_info.area_scr;
        loc_result.area_screening_enabled         = true;
        loc_result.area_screening_applied         = as.screen_applied;
        loc_result.area_screening_fallback_used   = as.fallback_used;
        loc_result.area_screening_fallback_reason = as.fallback_reason;
        loc_result.area_screening_center_xy       = as.center_xy_phys;
        loc_result.area_screening_radius_m        = as.radius_phys_m;
        loc_result.area_screening_pair_dist12_m   = as.pair_dist12_phys_m;
        loc_result.area_screening_candidate_xy    = as.pre_xy_phys;
        loc_result.area_screening_candidate_dist_to_center_m = as.dist_to_center_phys_m;
        loc_result.area_screening_pre_idx_shape   = as.pre_idx_shape;
        loc_result.area_screening_pre_dist_shape  = as.pre_dist_shape;
        loc_result.area_screening_valid_idx       = as.valid_idx;
        loc_result.area_screening_rejected_idx    = as.rejected_idx;
        loc_result.area_screening_valid_count     = as.valid_count;
    end
end


%% ==================== 局部函数 ====================

function true_pos = extract_true_position(event, FrameStates)
% EXTRACT_TRUE_POSITION  从 FrameStates 中提取事件对应的真值位置
%   取事件中间帧的真值位置

    true_pos = [];
    if isempty(FrameStates)
        return;
    end

    % 取事件中间帧
    t_mid = round((event.t_start + event.t_end) / 2);
    t_mid = max(1, min(t_mid, numel(FrameStates)));

    fs = FrameStates{t_mid};
    b  = event.band_id;

    if b <= numel(fs.active_per_band) && fs.active_per_band(b).has_source
        true_pos = fs.active_per_band(b).true_pos_xy;
    end
end
