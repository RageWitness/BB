function LocResults = run_m4_wknn_localization( ...
    EventList, SpatialFP, FrameStates, Config)
% RUN_M4_WKNN_LOCALIZATION  M4 模块总入口 — WKNN 定位
%
%   支持四种距离模式：
%     'shape_scale'         - 原始 shape+resid
%     'shape_scale_prob'    - 全局 AP 权重 + 概率 shape
%     'shape_scale_hn'      - 局部 hard-negative AP 权重 + 排斥项
%     'shape_scale_masked'  - 可靠 AP 掩码下的 shape+resid

    fprintf('\n============================================\n');
    fprintf('  M4 WKNN 定位模块\n');
    fprintf('============================================\n\n');

    tic;
    m4cfg = fill_m4_defaults(Config);
    fprintf('[M4] distance_mode = %s\n', m4cfg.distance_mode);
    fprintf('[M4] fingerprint_type = %s, fp_distance = %s\n', ...
        m4cfg.fingerprint_type, m4cfg.fp_distance);
    use_fp_l2 = ~strcmp(m4cfg.fingerprint_type, 'legacy');

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
        [F_obs_lin, F_obs_shape, F_obs_dBm] = aggregate_event_fingerprint_m4(ev);

        %% b. 全库 d_shape（用于 Area Screening 和诊断）
        F_lib_shape = SpatialFP.band(b).F_shape_l1;
        diff_shape  = F_obs_shape - F_lib_shape;
        d_shape_vec_full = sqrt(sum(diff_shape.^2, 1));   % 1 x G

        %% c. Area Screening
        scr = screen_candidates_by_area_m4( ...
            d_shape_vec_full, SpatialFP.grid_xy, m4cfg.area_screening);

        %% d. 构建 match_cfg
        match_cfg = struct();
        match_cfg.distance_mode   = m4cfg.distance_mode;
        match_cfg.lambda_shape    = m4cfg.lambda_shape;
        match_cfg.lambda_resid    = m4cfg.lambda_resid;
        match_cfg.lambda_hn       = m4cfg.lambda_hn;
        match_cfg.hn_delta_margin = m4cfg.hn_delta_margin;
        match_cfg.prob_eps        = m4cfg.prob_eps;
        match_cfg.masked_shape    = m4cfg.masked_shape;
        match_cfg.F_obs_dBm       = F_obs_dBm;

        [noise_floor_dBm, nf_ok, nf_msg] = get_noise_floor_dBm_for_band( ...
            Config, b, numel(F_obs_lin));
        if nf_ok
            match_cfg.noise_floor_dBm = noise_floor_dBm;
        else
            match_cfg.noise_floor_dBm = [];
            if strcmp(m4cfg.distance_mode, 'shape_scale_masked')
                fprintf('[M4] 警告: masked 模式噪声底缺失，自动回退原始距离 (%s)\n', nf_msg);
            end
        end

        %% e. 距离计算（分模式）
        if use_fp_l2
            % ---- 新框架：FP L1/L2 模式 ----
            preselect_idx = [];
            match_cfg_fp = struct();
            match_cfg_fp.fingerprint_type = m4cfg.fingerprint_type;
            match_cfg_fp.fp_distance      = m4cfg.fp_distance;

            G = numel(d_shape_vec_full);
            if scr.screen_applied
                vidx = scr.valid_idx;
                sub_fp = extract_band_subset_fp(SpatialFP.band(b), vidx);
                [D_sub, fp_extra] = compute_wknn_distance_fp_l2( ...
                    F_obs_lin, F_obs_shape, F_obs_dBm, sub_fp, match_cfg_fp);

                K_final = min(m4cfg.area_screening.k_final, numel(vidx));
                [nn_local, nn_dist, nn_w] = select_knn_neighbors_m4( ...
                    D_sub, K_final, m4cfg.eps_val);
                neighbor_idx  = vidx(nn_local);
                neighbor_dist = nn_dist;
                weights       = nn_w;

                d_shape_vec = d_shape_vec_full;
                q_opt_vec   = zeros(1, G);
                resid_vec   = inf(1, G);
                D_full      = inf(1, G);
                D_full(vidx) = D_sub;
            else
                [D_full, fp_extra] = compute_wknn_distance_fp_l2( ...
                    F_obs_lin, F_obs_shape, F_obs_dBm, SpatialFP.band(b), match_cfg_fp);
                [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4( ...
                    D_full, m4cfg.K, m4cfg.eps_val);
                d_shape_vec = d_shape_vec_full;
                q_opt_vec   = zeros(1, G);
                resid_vec   = inf(1, G);
            end

            prob_info = struct();
            prob_info.fp_extra = fp_extra;
            prob_info.D_vec_fp = D_full;

        elseif strcmp(m4cfg.distance_mode, 'shape_scale_hn')
            % ---- HN 模式：preselect → 子集距离 → WKNN ----
            [neighbor_idx, neighbor_dist, weights, ...
             d_shape_vec, q_opt_vec, resid_vec, prob_info, preselect_idx] = ...
                run_hn_pipeline(F_obs_lin, F_obs_shape, SpatialFP.band(b), ...
                    d_shape_vec_full, scr, m4cfg, match_cfg);
        else
            % ---- 原始 / prob / masked 模式 ----
            preselect_idx = [];
            if scr.screen_applied
                vidx = scr.valid_idx;
                sub_fp = extract_band_subset(SpatialFP.band(b), vidx, m4cfg.distance_mode);
                [D_sub, d_shape_sub, resid_sub, q_opt_sub, prob_info_sub] = ...
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
                d_shape_vec(vidx) = d_shape_sub;
                q_opt_vec   = zeros(1, G);  q_opt_vec(vidx) = q_opt_sub;
                resid_vec   = inf(1, G);    resid_vec(vidx) = resid_sub;
                prob_info   = prob_info_sub;

                if ~isempty(fieldnames(prob_info_sub)) && isfield(prob_info_sub, 'd_shape_prob_vec')
                    d_sp = inf(1, G); d_sp(vidx) = prob_info_sub.d_shape_prob_vec;
                    r_sp = inf(1, G); r_sp(vidx) = prob_info_sub.resid_prob_vec;
                    prob_info.d_shape_prob_vec = d_sp;
                    prob_info.resid_prob_vec   = r_sp;
                end
                if ~isempty(fieldnames(prob_info_sub)) && isfield(prob_info_sub, 'd_shape_mask_vec')
                    d_sm = inf(1, G); d_sm(vidx) = prob_info_sub.d_shape_mask_vec;
                    r_sm = inf(1, G); r_sm(vidx) = prob_info_sub.resid_mask_vec;
                    prob_info.d_shape_mask_vec = d_sm;
                    prob_info.resid_mask_vec   = r_sm;
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

        if isfield(prob_info, 'valid_ap_mask')
            extra_info.valid_ap_mask = prob_info.valid_ap_mask;
        end
        if isfield(prob_info, 'n_valid_ap')
            extra_info.n_valid_ap = prob_info.n_valid_ap;
        end

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
    [D_sub, ~, resid_sub, q_opt_sub, prob_info_sub] = ...
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
    resid_vec   = inf(1, G);    resid_vec(preselect_idx) = resid_sub;

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
        sub_fp.ap_weight_hn     = band_fp.ap_weight_hn(:, vidx);

        % hardneg_idx 存的是全库索引，不做子集映射，但需要
        % 指向全库的 mu/var，这里存全库引用供 compute_hn_mode 使用
        sub_fp.hardneg_idx = band_fp.hardneg_idx(vidx, :);
        sub_fp.full_F_shape_prob_mu  = band_fp.F_shape_prob_mu;
        sub_fp.full_F_shape_prob_var = band_fp.F_shape_prob_var;
    end
end


function sub_fp = extract_band_subset_fp(band_fp, vidx)
% EXTRACT_BAND_SUBSET_FP  为 FP L1/L2 模式提取候选子集
    sub_fp = struct();
    if isfield(band_fp, 'F_lin'),        sub_fp.F_lin        = band_fp.F_lin(:, vidx); end
    if isfield(band_fp, 'RF_raw'),       sub_fp.RF_raw       = band_fp.RF_raw(:, vidx); end
    if isfield(band_fp, 'RF_minmax'),    sub_fp.RF_minmax    = band_fp.RF_minmax(:, vidx); end
    if isfield(band_fp, 'F_shape_l1'),   sub_fp.F_shape_l1   = band_fp.F_shape_l1(:, vidx); end
    if isfield(band_fp, 'centered_dBm'), sub_fp.centered_dBm = band_fp.centered_dBm(:, vidx); end
    if isfield(band_fp, 'mean_dBm'),     sub_fp.mean_dBm     = band_fp.mean_dBm(vidx); end
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
        cfg.distance_mode = 'shape_scale_masked';
    end

    % 新框架：FP L1/L2 模式开关
    cfg.fingerprint_type = get_or_default(m4, 'fingerprint_type', 'legacy');
    cfg.fp_distance      = get_or_default(m4, 'fp_distance', 'L2');

    % K 近邻数
    cfg.K = get_or_default(m4, 'K', 3);
    cfg.eps_val = get_or_default(m4, 'eps_val', 1e-6);

    % lambda（根据模式设不同默认值）
    switch cfg.distance_mode
        case 'shape_scale_hn'
            def_ls = 0.50; def_lr = 0.25;
        case 'shape_scale_prob'
            def_ls = 0.60; def_lr = 0.40;
        case 'shape_scale_masked'
            def_ls = 0.70; def_lr = 0.30;
        otherwise
            def_ls = 0.70; def_lr = 0.30;
    end
    cfg.lambda_shape = get_or_default(m4, 'lambda_shape', def_ls);
    cfg.lambda_resid = get_or_default(m4, 'lambda_resid', def_lr);

    % prob 模式防除零
    if isfield(m4, 'prob_shape') && isfield(m4.prob_shape, 'eps')
        cfg.prob_eps = m4.prob_shape.eps;
    else
        cfg.prob_eps = 1e-12;
    end

    % masked 模式参数
    ms_defaults = struct();
    ms_defaults.enable = true;
    ms_defaults.tau_low_dB = 3;
    ms_defaults.tau_high_dB = 10;
    ms_defaults.eps = 1e-12;
    if isfield(m4, 'masked_shape')
        ms_in = m4.masked_shape;
        fns_ms = fieldnames(ms_defaults);
        for i = 1:numel(fns_ms)
            if isfield(ms_in, fns_ms{i})
                ms_defaults.(fns_ms{i}) = ms_in.(fns_ms{i});
            end
        end
    end
    cfg.masked_shape = ms_defaults;

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
    as_defaults.enable       = false;   % 保持旧逻辑默认
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


function [noise_floor_dBm, ok, reason] = get_noise_floor_dBm_for_band(Config, band_id, M)
% GET_NOISE_FLOOR_DBM_FOR_BAND  根据 n0 和 BW 返回当前频带噪声底 (dBm)
    if nargin < 3 || isempty(M)
        M = 1;
    end

    ok = false;
    reason = '';
    noise_floor_dBm = [];

    % 模式一：直接给每个 band 的噪声功率
    if isfield(Config, 'm1') && isfield(Config.m1, 'noise') && ...
       isfield(Config.m1.noise, 'mode') && strcmp(Config.m1.noise.mode, 'power_dBm')
        if isfield(Config.m1.noise, 'noise_power_dBm') && ...
           numel(Config.m1.noise.noise_power_dBm) >= band_id
            nf = Config.m1.noise.noise_power_dBm(band_id);
            noise_floor_dBm = repmat(nf, M, 1);
            ok = true;
            return;
        end
        reason = 'missing noise_power_dBm';
        return;
    end

    % 模式二：n0_dBmHz + BW
    if ~(isfield(Config, 'm1') && isfield(Config.m1, 'noise') && ...
         isfield(Config.m1.noise, 'n0_dBmHz'))
        reason = 'missing n0_dBmHz';
        return;
    end
    if ~(isfield(Config, 'm1') && isfield(Config.m1, 'bands') && ...
         isfield(Config.m1.bands, 'bw_Hz'))
        reason = 'missing bands.bw_Hz';
        return;
    end
    if numel(Config.m1.bands.bw_Hz) < band_id
        reason = 'bw_Hz index out of range';
        return;
    end

    n0_arr = Config.m1.noise.n0_dBmHz;
    if isscalar(n0_arr)
        n0 = n0_arr;
    elseif numel(n0_arr) >= band_id
        n0 = n0_arr(band_id);
    else
        reason = 'n0_dBmHz index out of range';
        return;
    end

    bw_hz = Config.m1.bands.bw_Hz(band_id);
    nf = n0 + 10 * log10(bw_hz);
    noise_floor_dBm = repmat(nf, M, 1);
    ok = true;
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
            errors(end+1) = lr.loc_error;
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
