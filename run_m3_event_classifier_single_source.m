function EventList = run_m3_event_classifier_single_source( ...
    Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config)
% RUN_M3_EVENT_CLASSIFIER_SINGLE_SOURCE  M3 模块总入口 — 单源事件分类
%
%   EventList = run_m3_event_classifier_single_source(
%       Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config)
%
%   流程：
%     1. detect_band_events_m3      — 逐频带事件检测
%     2. extract_event_features_m3  — 事件特征提取
%     3. score_event_type_m3        — 四类源评分与判定
%     4. route_event_m3             — 路由分配
%
%   输入：
%       Y_dBm_all       - (M x B x T) 全时段 RSS 观测 (dBm)
%       Y_lin_all       - (M x B x T) 全时段 RSS 观测 (线性 mW)
%       SpatialFP       - 空间指纹库
%       SignatureLib    - 静态模板特征库
%       SourceTemplates - M0 源模板
%       Config          - 配置
%
%   输出：
%       EventList(e) — 完整事件列表，每个元素包含：
%         基本信息：event_id, band_id, time_range, t_start, t_end, duration
%         观测数据：obs_segment_dBm (MxL), obs_segment_lin (MxL)
%         特征：feat_vec, power_level_est, power_stability_est,
%               band_coverage_vec, schedule_match_score, position_prior_match_score
%         评分：score_trusted, score_prior_pos, score_prior_time, score_target
%         判定：type_hat, route_action, linked_template_key, upgrade_hint

    fprintf('\n============================================\n');
    fprintf('  M3 单源事件分类模块\n');
    fprintf('============================================\n\n');

    tic;

    %% 1. 事件检测
    EventListRaw = detect_band_events_m3(Y_dBm_all, Y_lin_all, Config);

    if isempty(EventListRaw)
        fprintf('[M3] 未检测到任何事件\n');
        EventList = [];
        return;
    end

    %% 2. 事件特征提取
    EventList = extract_event_features_m3(EventListRaw, Y_dBm_all, Y_lin_all, ...
        SpatialFP, SignatureLib, SourceTemplates, Config);

    %% 3. 四类源评分与判定
    EventList = score_event_type_m3(EventList, Config);

    %% 4. 路由分配
    EventList = route_event_m3(EventList, Config);

    elapsed = toc;
    fprintf('\n[M3] 完成: %d 个事件, 耗时 %.2f s\n', numel(EventList), elapsed);

    %% 打印事件摘要
    print_event_summary(EventList);

    fprintf('\n============================================\n');
    fprintf('  M3 事件分类完成\n');
    fprintf('============================================\n');
end


%% ==================== 局部函数 ====================

function print_event_summary(EventList)
% PRINT_EVENT_SUMMARY  打印事件分类摘要表

    n = numel(EventList);
    if n == 0, return; end

    fprintf('\n--- 事件分类摘要 (前 %d 个) ---\n', min(n, 20));
    fprintf('  %-6s %-6s %-10s %-6s %-8s %-12s %-25s %-12s\n', ...
        'ID', 'Band', 'Time', 'Dur', 'Power', 'Type', 'Scores[tr/pp/pt/tg]', 'Route');

    for e = 1:min(n, 20)
        ev = EventList(e);
        fprintf('  %-6d %-6d [%3d,%3d]  %-6d %-8.1f %-12s [%.2f/%.2f/%.2f/%.2f]  %-12s\n', ...
            ev.event_id, ev.band_id, ev.t_start, ev.t_end, ev.duration, ...
            ev.power_level_est, ev.type_hat, ...
            ev.score_trusted, ev.score_prior_pos, ev.score_prior_time, ev.score_target, ...
            ev.route_action);
    end
    if n > 20
        fprintf('  ... 共 %d 个事件 (仅显示前 20 个)\n', n);
    end
end
