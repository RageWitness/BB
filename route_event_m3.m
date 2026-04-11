function EventList = route_event_m3(EventList, Config)
% ROUTE_EVENT_M3  根据事件类型判定结果决定路由动作
%
%   EventList = route_event_m3(EventList, Config)
%
%   路由规则：
%     trusted_fixed     → 'calibrate_direct'
%     prior_pos_known   → 'calibrate_direct'
%     prior_time_known  → 'localize_then_calibrate'
%     ordinary_target   → 'localize_only'
%     若最高分太低或与次高分差距过小 → 'hold'
%
%   同时输出：
%     linked_template_key — 最佳匹配模板 key（对 trusted/prior 类有意义）
%     upgrade_hint        — 对 ordinary_target 的晋升提示
%
%   输入/输出：
%       EventList(e) — 增加 route_action, linked_template_key, upgrade_hint

    % --- 填充路由参数 ---
    rt = fill_route_defaults(Config);

    n_events = numel(EventList);

    % 预分配新字段
    for e = 1:n_events
        EventList(e).route_action        = '';
        EventList(e).linked_template_key = '';
        EventList(e).upgrade_hint        = '';
    end

    route_counts = struct('calibrate_direct', 0, ...
                          'localize_then_calibrate', 0, ...
                          'localize_only', 0, ...
                          'hold', 0);

    for e = 1:n_events
        ev = EventList(e);

        scores = [ev.score_trusted, ev.score_prior_pos, ...
                  ev.score_prior_time, ev.score_target];

        [max_score, best_idx] = max(scores);
        sorted_scores = sort(scores, 'descend');
        margin = sorted_scores(1) - sorted_scores(2);

        %% 判断是否进入 hold
        if max_score < rt.min_score_threshold || margin < rt.min_margin
            ev.route_action = 'hold';
            ev.linked_template_key = '';
            ev.upgrade_hint = '';
            route_counts.hold = route_counts.hold + 1;
            EventList(e) = ev;
            continue;
        end

        %% 根据 type_hat 设定路由
        switch ev.type_hat
            case 'trusted_fixed'
                ev.route_action = 'calibrate_direct';
                route_counts.calibrate_direct = route_counts.calibrate_direct + 1;

            case 'prior_pos_known'
                ev.route_action = 'calibrate_direct';
                route_counts.calibrate_direct = route_counts.calibrate_direct + 1;

            case 'prior_time_known'
                ev.route_action = 'localize_then_calibrate';
                route_counts.localize_then_calibrate = route_counts.localize_then_calibrate + 1;

            case 'ordinary_target'
                ev.route_action = 'localize_only';
                route_counts.localize_only = route_counts.localize_only + 1;

            otherwise
                ev.route_action = 'hold';
                route_counts.hold = route_counts.hold + 1;
        end

        %% linked_template_key（暂置空，后续可精确关联）
        ev.linked_template_key = '';

        %% upgrade_hint
        if strcmp(ev.type_hat, 'ordinary_target')
            % 根据功率稳定性和持续时间给出晋升提示
            if ev.power_stability_est < rt.upgrade_stability_thresh && ...
               ev.duration >= rt.upgrade_duration_thresh
                ev.upgrade_hint = 'high';
            elseif ev.power_stability_est < rt.upgrade_stability_thresh * 2
                ev.upgrade_hint = 'medium';
            else
                ev.upgrade_hint = 'low';
            end
        else
            ev.upgrade_hint = 'none';
        end

        EventList(e) = ev;
    end

    fprintf('[M3] 路由分配: calibrate_direct=%d, loc_then_cal=%d, loc_only=%d, hold=%d\n', ...
        route_counts.calibrate_direct, ...
        route_counts.localize_then_calibrate, ...
        route_counts.localize_only, ...
        route_counts.hold);
end


%% ==================== 局部函数 ====================

function rt = fill_route_defaults(Config)
% FILL_ROUTE_DEFAULTS  填充路由参数默认值

    if isfield(Config, 'm3') && isfield(Config.m3, 'route')
        r = Config.m3.route;
    else
        r = struct();
    end

    % 最低分数阈值：低于此值的事件进入 hold
    if isfield(r, 'min_score_threshold')
        rt.min_score_threshold = r.min_score_threshold;
    else
        rt.min_score_threshold = 0.15;
    end

    % 最高与次高分差距阈值
    if isfield(r, 'min_margin')
        rt.min_margin = r.min_margin;
    else
        rt.min_margin = 0.05;
    end

    % 晋升提示阈值
    if isfield(r, 'upgrade_stability_thresh')
        rt.upgrade_stability_thresh = r.upgrade_stability_thresh;
    else
        rt.upgrade_stability_thresh = 2.0;  % dB^2
    end

    if isfield(r, 'upgrade_duration_thresh')
        rt.upgrade_duration_thresh = r.upgrade_duration_thresh;
    else
        rt.upgrade_duration_thresh = 5;  % 帧
    end
end
