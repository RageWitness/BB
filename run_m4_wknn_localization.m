function LocResults = run_m4_wknn_localization( ...
    EventList, SpatialFP, FrameStates, Config)
% RUN_M4_WKNN_LOCALIZATION  M4 模块总入口 — 简易 WKNN 定位
%
%   LocResults = run_m4_wknn_localization(
%       EventList, SpatialFP, FrameStates, Config)
%
%   只对需要定位的事件执行 WKNN：
%     route_action = 'localize_only'           → ordinary_target
%     route_action = 'localize_then_calibrate' → prior_time_known
%
%   流程：
%     1. 筛选待定位事件
%     2. 对每个事件：
%        a. aggregate_event_fingerprint_m4     — 聚合指纹
%        b. compute_wknn_distance_power_corrected — 功率修正距离
%        c. select_knn_neighbors_m4            — 选 K 近邻
%        d. estimate_position_wknn_m4          — 加权位置估计
%        e. record_loc_results_m4              — 记录结果
%
%   输入：
%       EventList   - M3 输出的事件列表
%       SpatialFP   - 空间指纹库
%       FrameStates - cell(1,T) 帧状态（含真值）
%       Config      - 配置
%
%   输出：
%       LocResults(k) — 定位结果列表

    fprintf('\n============================================\n');
    fprintf('  M4 WKNN 定位模块\n');
    fprintf('============================================\n\n');

    tic;

    % --- 填充默认参数 ---
    m4cfg = fill_m4_defaults(Config);

    % --- 筛选待定位事件 ---
    if isempty(EventList)
        fprintf('[M4] 无事件输入\n');
        LocResults = [];
        return;
    end

    loc_mask = false(1, numel(EventList));
    for e = 1:numel(EventList)
        ra = EventList(e).route_action;
        if strcmp(ra, 'localize_only') || strcmp(ra, 'localize_then_calibrate')
            loc_mask(e) = true;
        end
    end
    loc_events = EventList(loc_mask);
    n_loc = numel(loc_events);

    fprintf('[M4] 待定位事件: %d / %d\n', n_loc, numel(EventList));

    if n_loc == 0
        LocResults = [];
        return;
    end

    LocResults = [];

    for k = 1:n_loc
        ev = loc_events(k);
        b  = ev.band_id;

        %% a. 聚合事件指纹（含形状向量）
        [F_obs_lin, F_obs_shape] = aggregate_event_fingerprint_m4(ev);

        %% b. 全库形状距离（仅 d_shape，用于 Area Screening 预候选）
        match_cfg = struct();
        match_cfg.lambda_shape = m4cfg.lambda_shape;
        match_cfg.lambda_resid = m4cfg.lambda_resid;

        F_lib_shape = SpatialFP.band(b).F_shape_l1;            % M x G
        diff_shape  = F_obs_shape - F_lib_shape;                % M x G
        d_shape_vec_full = sqrt(sum(diff_shape.^2, 1));         % 1 x G

        %% c. Area Screening — 空间一致性筛选
        scr = screen_candidates_by_area_m4( ...
            d_shape_vec_full, SpatialFP.grid_xy, m4cfg.area_screening);

        if scr.screen_applied
            %% —— 筛选生效：仅对 valid 候选计算 q_g*, resid, D ——
            vidx = scr.valid_idx;                               % 通过筛选的网格索引

            % 构造子集指纹库
            sub_fp = struct();
            sub_fp.F_lin      = SpatialFP.band(b).F_lin(:, vidx);
            sub_fp.F_shape_l1 = SpatialFP.band(b).F_shape_l1(:, vidx);

            [D_sub, d_shape_sub, resid_sub, q_opt_sub] = ...
                compute_wknn_distance_shape_scale( ...
                    F_obs_lin, F_obs_shape, sub_fp, match_cfg);

            % 最终 K 近邻数
            K_final = min(m4cfg.area_screening.k_final, numel(vidx));

            [nn_local_idx, nn_dist, nn_weights] = select_knn_neighbors_m4( ...
                D_sub, K_final, m4cfg.eps_val);

            % 映射回全局索引
            neighbor_idx  = vidx(nn_local_idx);
            neighbor_dist = nn_dist;
            weights       = nn_weights;

            % 还原全库向量（用于诊断记录）
            G = numel(d_shape_vec_full);
            d_shape_vec = d_shape_vec_full;
            q_opt_vec   = zeros(1, G);  q_opt_vec(vidx)   = q_opt_sub;
            resid_vec   = inf(1, G);    resid_vec(vidx)    = resid_sub;
            D_vec       = inf(1, G);    D_vec(vidx)        = D_sub;

        else
            %% —— 回退到原始全库逻辑 ——
            [D_vec, d_shape_vec, resid_vec, q_opt_vec] = ...
                compute_wknn_distance_shape_scale( ...
                    F_obs_lin, F_obs_shape, SpatialFP.band(b), match_cfg);

            [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4( ...
                D_vec, m4cfg.K, m4cfg.eps_val);
        end

        %% d. 加权位置估计
        est_pos_xy = estimate_position_wknn_m4( ...
            SpatialFP.grid_xy, neighbor_idx, weights);

        %% e. 记录结果（含诊断信息 + Area Screening 诊断）
        extra_info = struct();
        extra_info.q_opt_vec   = q_opt_vec;
        extra_info.d_shape_vec = d_shape_vec;
        extra_info.resid_vec   = resid_vec;
        extra_info.area_scr    = scr;

        lr = record_loc_results_m4(ev, est_pos_xy, F_obs_lin, ...
            neighbor_idx, neighbor_dist, weights, FrameStates, Config, extra_info);

        if isempty(LocResults)
            LocResults = lr;
        else
            LocResults(end+1) = lr; %#ok<AGROW>
        end
    end

    elapsed = toc;
    fprintf('\n[M4] 完成: %d 个事件定位, 耗时 %.2f s\n', n_loc, elapsed);

    %% 打印定位结果摘要
    print_loc_summary(LocResults);

    fprintf('\n============================================\n');
    fprintf('  M4 WKNN 定位完成\n');
    fprintf('============================================\n');
end


%% ==================== 局部函数 ====================

function cfg = fill_m4_defaults(Config)
% FILL_M4_DEFAULTS  填充 M4 默认参数（含 Area Screening）

    if isfield(Config, 'm4')
        m4 = Config.m4;
    else
        m4 = struct();
    end

    % K 近邻数
    if isfield(m4, 'K')
        cfg.K = m4.K;
    else
        cfg.K = 3;
    end

    % 防除零常数
    if isfield(m4, 'eps_val')
        cfg.eps_val = m4.eps_val;
    else
        cfg.eps_val = 1e-6;
    end

    % 形状+缩放距离权重
    if isfield(m4, 'lambda_shape')
        cfg.lambda_shape = m4.lambda_shape;
    else
        cfg.lambda_shape = 0.7;
    end
    if isfield(m4, 'lambda_resid')
        cfg.lambda_resid = m4.lambda_resid;
    else
        cfg.lambda_resid = 0.3;
    end

    % ---- Area Screening 默认参数 ----
    as_defaults = struct();
    as_defaults.enable       = true;
    as_defaults.k_shape      = 5;
    as_defaults.k_final      = 6;
    as_defaults.r_min_m      = 8;
    as_defaults.alpha        = 1.5;
    as_defaults.r_max_m      = 15;
    as_defaults.d_pair_max_m = 40;
    as_defaults.k_min_valid  = 2;

    if isfield(m4, 'area_screening')
        as_in = m4.area_screening;
        fns = fieldnames(as_defaults);
        for i = 1:numel(fns)
            if isfield(as_in, fns{i})
                as_defaults.(fns{i}) = as_in.(fns{i});
            end
        end
    end
    cfg.area_screening = as_defaults;
end


function print_loc_summary(LocResults)
% PRINT_LOC_SUMMARY  打印定位结果摘要

    n = numel(LocResults);
    if n == 0, return; end

    fprintf('\n--- 定位结果摘要 ---\n');
    fprintf('  %-6s %-6s %-18s %-18s %-12s %-10s\n', ...
        'ID', 'Band', 'Type', 'Est_pos', 'Error(m)', 'Route');

    errors = [];
    for k = 1:n
        lr = LocResults(k);
        if ~isnan(lr.loc_error)
            errors(end+1) = lr.loc_error; %#ok<AGROW>
        end
        fprintf('  %-6d %-6d %-18s [%6.1f,%6.1f]  %-12s %-10s\n', ...
            lr.event_id, lr.band_id, lr.type_hat, ...
            lr.est_pos_xy(1), lr.est_pos_xy(2), ...
            format_error(lr.loc_error), lr.route_action);
    end

    if ~isempty(errors)
        fprintf('\n--- 定位误差统计 ---\n');
        fprintf('  事件数: %d (有真值: %d)\n', n, numel(errors));
        fprintf('  平均误差: %.2f m\n', mean(errors));
        fprintf('  中位误差: %.2f m\n', median(errors));
        fprintf('  RMSE:     %.2f m\n', sqrt(mean(errors.^2)));
        fprintf('  最大误差: %.2f m\n', max(errors));
        fprintf('  最小误差: %.2f m\n', min(errors));

        % 按类型分别统计
        types_loc = {LocResults.type_hat};
        for type_cell = {'ordinary_target', 'prior_time_known'}
            type_name = type_cell{1};
            mask = strcmp(types_loc, type_name);
            errs = [LocResults(mask).loc_error];
            errs = errs(~isnan(errs));
            if ~isempty(errs)
                fprintf('  [%s] n=%d, mean=%.2f m, RMSE=%.2f m\n', ...
                    type_name, numel(errs), mean(errs), sqrt(mean(errs.^2)));
            end
        end
    end
end


function s = format_error(err)
% FORMAT_ERROR  格式化误差输出
    if isnan(err)
        s = 'N/A';
    else
        s = sprintf('%.2f', err);
    end
end
