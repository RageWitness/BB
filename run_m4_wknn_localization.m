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

        %% a. 聚合事件指纹
        F_obs = aggregate_event_fingerprint_m4(ev);

        %% b. 功率修正距离计算
        [dist_vec, ~] = compute_wknn_distance_power_corrected( ...
            F_obs, SpatialFP.band(b));

        %% c. 选 K 近邻
        [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4( ...
            dist_vec, m4cfg.K, m4cfg.eps_val);

        %% d. 加权位置估计
        est_pos_xy = estimate_position_wknn_m4( ...
            SpatialFP.grid_xy, neighbor_idx, weights);

        %% e. 记录结果
        lr = record_loc_results_m4(ev, est_pos_xy, F_obs, ...
            neighbor_idx, neighbor_dist, weights, FrameStates, Config);

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
% FILL_M4_DEFAULTS  填充 M4 默认参数

    if isfield(Config, 'm4')
        m4 = Config.m4;
    else
        m4 = struct();
    end

    % K 近邻数
    if isfield(m4, 'K')
        cfg.K = m4.K;
    else
        cfg.K = 5;
    end

    % 防除零常数
    if isfield(m4, 'eps_val')
        cfg.eps_val = m4.eps_val;
    else
        cfg.eps_val = 1e-6;
    end
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
