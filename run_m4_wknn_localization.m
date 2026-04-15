function LocResults = run_m4_wknn_localization( ...
    EventList, SpatialFP, FrameStates, Config)
% RUN_M4_WKNN_LOCALIZATION  M4 模块总入口 — WKNN 定位
%
%   支持三种距离模式：
%     'shape_scale'       — 原始 shape+resid
%     'shape_scale_prob'  — 全局 AP 权重 + 概率 shape
%     'shape_scale_hn'    — 局部 hard-negative AP 权重 + 排斥项

    fprintf('\n============================================\n');
    fprintf('  M4 WKNN 定位模块\n');
    fprintf('============================================\n\n');

    tic;
    m4cfg = fill_m4_defaults(Config);
    fprintf('[M4] distance_mode = %s\n', m4cfg.distance_mode);

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

        %% a. 聚合事件指纹
        [F_obs_lin, F_obs_shape] = aggregate_event_fingerprint_m4(ev);

        %% b. 全库 d_shape（用于 Area Screening 和预筛）
        F_lib_shape = SpatialFP.band(b).F_shape_l1;
        diff_shape  = F_obs_shape - F_lib_shape;
        d_shape_vec_full = sqrt(sum(diff_shape.^2, 1));   % 1 x G

        %% c. Area Screening
        scr = screen_candidates_by_area_m4( ...
            d_shape_vec_full, SpatialFP.grid_xy, m4cfg.area_screening);

        %% d. 构建 match_cfg
        match_cfg = struct();
        match_cfg.distance_mode  = m4cfg.distance_mode;
        match_cfg.lambda_shape   = m4cfg.lambda_shape;
        match_cfg.lambda_resid   = m4cfg.lambda_resid;
        match_cfg.lambda_hn      = m4cfg.lambda_hn;
        match_cfg.hn_delta_margin = m4cfg.hn_delta_margin;
        match_cfg.prob_eps       = m4cfg.prob_eps;

        %% e. 距离计算（分模式）
        if strcmp(m4cfg.distance_mode, 'shape_scale_hn')
            % ---- HN 模式：preselect → 子集距离 → WKNN ----
            [neighbor_idx, neighbor_dist, weights, ...
             d_shape_vec, q_opt_vec, resid_vec, prob_info, preselect_idx] = ...
                run_hn_pipeline(F_obs_lin, F_obs_shape, SpatialFP.band(b), ...
                    d_shape_vec_full, scr, m4cfg, match_cfg);
        else
            % ---- 原始 / prob 模式 ----
            preselect_idx = [];
            if scr.screen_applied
                vidx = scr.valid_idx;
                sub_fp = extract_band_subset(SpatialFP.band(b), vidx, m4cfg.distance_mode);
                [D_sub, ~, resid_sub, q_opt_sub, prob_info_sub] = ...
                    compute_wknn_distance_shape_scale( ...
                        F_obs_lin, F_obs_shape, sub_fp, match_cfg);
                K_final = min(m4cfg.area_screening.k_final, numel(vidx));
                [nn_local, nn_dist, nn_w] = select_knn_neighbors_m4( ...
                    D_sub, K_final, m4cfg.eps_val);
                neighbor_idx  = vidx(nn_local);
                neighbor_dist = nn_dist;
                weights       = nn_w;
                G = numel(d_shape_vec_full);
                d_shape_vec = d_shape_vec_full;
                q_opt_vec   = zeros(1, G);  q_opt_vec(vidx) = q_opt_sub;
                resid_vec   = inf(1, G);    resid_vec(vidx)  = resid_sub;
                prob_info   = prob_info_sub;
                if ~isempty(fieldnames(prob_info_sub)) && isfield(prob_info_sub, 'd_shape_prob_vec')
                    d_sp = inf(1, G);  d_sp(vidx) = prob_info_sub.d_shape_prob_vec;
                    r_sp = inf(1, G);  r_sp(vidx) = prob_info_sub.resid_prob_vec;
                    prob_info.d_shape_prob_vec = d_sp;
                    prob_info.resid_prob_vec   = r_sp;
                end
            else
                [D_vec, d_shape_vec, resid_vec, q_opt_vec, prob_info] = ...
                    compute_wknn_distance_shape_scale( ...
                        F_obs_lin, F_obs_shape, SpatialFP.band(b), match_cfg);
                [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4( ...
                    D_vec, m4cfg.K, m4cfg.eps_val);
            end
        end

        %% f. 加权位置估计
        est_pos_xy = estimate_position_wknn_m4( ...
            SpatialFP.grid_xy, neighbor_idx, weights);

        %% g. 记录结果
        extra_info = struct();
        extra_info.q_opt_vec     = q_opt_vec;
        extra_info.d_shape_vec   = d_shape_vec;
        extra_info.resid_vec     = resid_vec;
        extra_info.area_scr      = scr;
        extra_info.distance_mode = m4cfg.distance_mode;
        extra_info.prob_info     = prob_info;
        extra_info.preselect_idx = preselect_idx;

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
    print_loc_summary(LocResults);
    fprintf('\n============================================\n');
    fprintf('  M4 WKNN 定位完成\n');
    fprintf('============================================\n');
end


%% ==================== 局部函数 ====================

function [neighbor_idx, neighbor_dist, weights, ...
          d_shape_vec, q_opt_vec, resid_vec, prob_info, preselect_idx] = ...
    run_hn_pipeline(F_obs_lin, F_obs_shape, band_fp, d_shape_full, scr, m4cfg, match_cfg)
% RUN_HN_PIPELINE  shape_scale_hn 模式的完整管线
%
%   1. preselect：用原始 d_shape 取前 preselect_k 个候选
%      （若 area screening 生效，先从 valid 集合内选）
%   2. 对 preselect 子集计算 D_hn
%   3. 在子集上做最终 WKNN

    G = numel(d_shape_full);
    preselect_k = m4cfg.hn_preselect_k;

    % ---- 确定 preselect 候选池 ----
    if scr.screen_applied
        pool_idx = scr.valid_idx;
        d_pool   = d_shape_full(pool_idx);
    else
        pool_idx = 1:G;
        d_pool   = d_shape_full;
    end

    n_pool = numel(pool_idx);
    n_pre  = min(preselect_k, n_pool);
    [~, sort_order] = sort(d_pool, 'ascend');
    pre_local  = sort_order(1:n_pre);
    preselect_idx = pool_idx(pre_local);          % 全局索引

    % ---- 构造 preselect 子集指纹库 ----
    sub_fp = extract_band_subset(band_fp, preselect_idx, 'shape_scale_hn');

    % ---- 计算 D_hn ----
    [D_sub, d_shape_sub, resid_sub, q_opt_sub, prob_info_sub] = ...
        compute_wknn_distance_shape_scale( ...
            F_obs_lin, F_obs_shape, sub_fp, match_cfg);

    % ---- 最终 WKNN ----
    K_final = min(m4cfg.K, n_pre);
    [nn_local, nn_dist, nn_w] = select_knn_neighbors_m4( ...
        D_sub, K_final, m4cfg.eps_val);

    neighbor_idx  = preselect_idx(nn_local);
    neighbor_dist = nn_dist;
    weights       = nn_w;

    % ---- 还原全库向量 ----
    d_shape_vec = d_shape_full;
    q_opt_vec   = zeros(1, G);  q_opt_vec(preselect_idx) = q_opt_sub;
    resid_vec   = inf(1, G);    resid_vec(preselect_idx)  = resid_sub;

    % 还原 hn_info 向量到全库大小
    prob_info = prob_info_sub;
    if isfield(prob_info_sub, 'd_shape_hn_vec')
        flds = {'d_shape_hn_vec', 'resid_prob_vec', 'd_neg_min_vec', 'p_hn_vec'};
        for fi = 1:numel(flds)
            fn = flds{fi};
            full_v = inf(1, G);
            if strcmp(fn, 'p_hn_vec')
                full_v = zeros(1, G);
            end
            full_v(preselect_idx) = prob_info_sub.(fn);
            prob_info.(fn) = full_v;
        end
    end
end


function sub_fp = extract_band_subset(band_fp, vidx, distance_mode)
% EXTRACT_BAND_SUBSET  从频带指纹库中提取子集

    sub_fp = struct();
    sub_fp.F_lin      = band_fp.F_lin(:, vidx);
    sub_fp.F_shape_l1 = band_fp.F_shape_l1(:, vidx);

    if strcmp(distance_mode, 'shape_scale_prob')
        sub_fp.F_shape_prob_mu  = band_fp.F_shape_prob_mu(:, vidx);
        sub_fp.F_shape_prob_var = band_fp.F_shape_prob_var(:, vidx);
        sub_fp.ap_weight_global = band_fp.ap_weight_global;
        sub_fp.std_lin          = band_fp.std_lin(:, vidx);

    elseif strcmp(distance_mode, 'shape_scale_hn')
        sub_fp.F_shape_prob_mu  = band_fp.F_shape_prob_mu(:, vidx);
        sub_fp.F_shape_prob_var = band_fp.F_shape_prob_var(:, vidx);
        sub_fp.std_lin          = band_fp.std_lin(:, vidx);
        sub_fp.ap_weight_hn    = band_fp.ap_weight_hn(:, vidx);

        % hardneg_idx 存的是全库索引，不做子集映射，但需要
        % 指向全库的 mu/var，这里存全库引用供 compute_hn_mode 使用
        sub_fp.hardneg_idx     = band_fp.hardneg_idx(vidx, :);

        % HN 模式需要全库 mu/var 来算 d_neg（因为 hardneg 指向全库点）
        sub_fp.full_F_shape_prob_mu  = band_fp.F_shape_prob_mu;
        sub_fp.full_F_shape_prob_var = band_fp.F_shape_prob_var;
    end
end


function cfg = fill_m4_defaults(Config)
% FILL_M4_DEFAULTS  填充 M4 默认参数

    if isfield(Config, 'm4')
        m4 = Config.m4;
    else
        m4 = struct();
    end

    % 距离模式
    if isfield(m4, 'distance_mode')
        cfg.distance_mode = m4.distance_mode;
    else
        cfg.distance_mode = 'shape_scale_hn';
    end

    % K 近邻数
    cfg.K = get_or_default(m4, 'K', 3);
    cfg.eps_val = get_or_default(m4, 'eps_val', 1e-6);

    % lambda（根据模式设不同默认值）
    switch cfg.distance_mode
        case 'shape_scale_hn'
            def_ls = 0.50; def_lr = 0.25;
        case 'shape_scale_prob'
            def_ls = 0.6;  def_lr = 0.4;
        otherwise
            def_ls = 0.7;  def_lr = 0.3;
    end
    cfg.lambda_shape = get_or_default(m4, 'lambda_shape', def_ls);
    cfg.lambda_resid = get_or_default(m4, 'lambda_resid', def_lr);

    % prob 模式防除零
    if isfield(m4, 'prob_shape') && isfield(m4.prob_shape, 'eps')
        cfg.prob_eps = m4.prob_shape.eps;
    else
        cfg.prob_eps = 1e-12;
    end

    % ---- HN 模式参数 ----
    if isfield(m4, 'hardneg')
        hn = m4.hardneg;
    else
        hn = struct();
    end
    cfg.hn_preselect_k = get_or_default(hn, 'preselect_k', 20);
    cfg.hn_delta_margin = get_or_default(hn, 'delta_margin', 0.15);
    cfg.lambda_hn      = get_or_default(hn, 'lambda_hn', 0.25);

    % ---- Area Screening 默认参数 ----
    as_defaults = struct();
    as_defaults.enable       = false;   % HN 模式默认关闭 area screening
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


function val = get_or_default(s, fn, def)
    if isfield(s, fn)
        val = s.(fn);
    else
        val = def;
    end
end


function print_loc_summary(LocResults)
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
        fprintf('  P95 误差: %.2f m\n', prctile(errors, 95));

        types_loc = {LocResults.type_hat};
        for type_cell = {'ordinary_target', 'prior_time_known'}
            type_name = type_cell{1};
            mask = strcmp(types_loc, type_name);
            errs = [LocResults(mask).loc_error];
            errs = errs(~isnan(errs));
            if ~isempty(errs)
                fprintf('  [%s] n=%d, mean=%.2f m, RMSE=%.2f m, P95=%.2f m\n', ...
                    type_name, numel(errs), mean(errs), sqrt(mean(errs.^2)), prctile(errs,95));
            end
        end
    end
end


function s = format_error(err)
    if isnan(err)
        s = 'N/A';
    else
        s = sprintf('%.2f', err);
    end
end
