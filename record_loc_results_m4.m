function loc_result = record_loc_results_m4(event, est_pos_xy, F_obs, ...
    neighbor_idx, neighbor_dist, weights, FrameStates, Config, extra_info)
% RECORD_LOC_RESULTS_M4  记录单个事件的定位结果
%
%   loc_result = record_loc_results_m4(event, est_pos_xy, F_obs,
%       neighbor_idx, neighbor_dist, weights, FrameStates, Config)
%   loc_result = record_loc_results_m4(event, est_pos_xy, F_obs,
%       neighbor_idx, neighbor_dist, weights, FrameStates, Config, extra_info)
%
%   输出 loc_result 包含：
%     event_id, band_id, type_hat, route_action
%     obs_fp_lin        — 事件聚合指纹 (M x 1, 线性 mW)
%     obs_fp_dBm        — 事件聚合指纹 (M x 1, dBm)
%     est_pos_xy        — 估计位置 (1 x 2)
%     neighbor_idx      — K 近邻索引
%     neighbor_dist     — K 近邻距离
%     weights           — 归一化权重
%     match_score       — 最佳匹配距离的倒数
%     true_pos_xy       — 真值位置（若可从 FrameStates 获取）
%     loc_error         — 定位误差 (m)
%     n0                — 当前噪声底 (dBm/Hz)
%     q_opt_neighbors   — (1xK) 近邻的最优缩放系数（若提供 extra_info）
%     d_shape_neighbors — (1xK) 近邻的形状距离（若提供 extra_info）
%     resid_neighbors   — (1xK) 近邻的重构残差（若提供 extra_info）
%
%   输入：
%       event           - 事件结构体
%       est_pos_xy      - (1 x 2) 估计位置
%       F_obs           - (M x 1) 聚合指纹
%       neighbor_idx    - (1 x K)
%       neighbor_dist   - (1 x K)
%       weights         - (1 x K)
%       FrameStates     - cell(1,T) 帧状态
%       Config          - 配置
%       extra_info      - [可选] 含 q_opt_vec, d_shape_vec, resid_vec

    loc_result = struct();

    % 基本信息
    loc_result.event_id     = event.event_id;
    loc_result.band_id      = event.band_id;
    loc_result.type_hat     = event.type_hat;
    loc_result.route_action = event.route_action;
    loc_result.time_range   = event.time_range;

    % 观测与估计
    loc_result.obs_fp_lin   = F_obs;
    loc_result.obs_fp_dBm   = 10 * log10(max(F_obs, 1e-30));
    loc_result.est_pos_xy   = est_pos_xy;

    % 近邻信息
    loc_result.neighbor_idx  = neighbor_idx;
    loc_result.neighbor_dist = neighbor_dist;
    loc_result.weights       = weights;

    % 匹配分数：最小距离的倒数（越大越好）
    loc_result.match_score = 1 / (min(neighbor_dist) + 1e-10);

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
        loc_result.n0 = -174;  % 默认热噪声底
    end

    % 形状+缩放诊断信息
    if nargin >= 9 && ~isempty(extra_info)
        loc_result.q_opt_neighbors   = extra_info.q_opt_vec(neighbor_idx);
        loc_result.d_shape_neighbors = extra_info.d_shape_vec(neighbor_idx);
        loc_result.resid_neighbors   = extra_info.resid_vec(neighbor_idx);
    else
        loc_result.q_opt_neighbors   = [];
        loc_result.d_shape_neighbors = [];
        loc_result.resid_neighbors   = [];
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
