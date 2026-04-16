function [template_scores, best_template_key, best_score, second_best_score, diag_info] = ...
    compute_group_time_match_m3(GroupItem, TemplateList, Config)
% COMPUTE_GROUP_TIME_MATCH_M3  group 级逐模板时间匹配
%
%   S_time_j = a1*IoU + a2*exp(-delta_start^2/(2*sigma_t^2)) ...
%            + a3*exp(-delta_duration^2/(2*sigma_d^2))
%   若周期模板，再乘以 phase 项。

    cfg = fill_time_match_defaults(Config);
    N = numel(TemplateList);

    template_scores = zeros(1, N);
    best_template_key = '';
    best_score = 0;
    second_best_score = 0;

    diag_info = struct();
    diag_info.iou_best = zeros(1, N);
    diag_info.delta_start_best = inf(1, N);
    diag_info.delta_duration_best = inf(1, N);
    diag_info.phase_score = zeros(1, N);

    if N == 0
        return;
    end

    obs_start = GroupItem.start_frame;
    obs_end = GroupItem.end_frame;
    obs_dur = max(1, obs_end - obs_start + 1);

    for j = 1:N
        tpl = TemplateList(j);
        tp = resolve_time_prior_local(tpl);
        [intervals, phase_delta] = build_template_intervals(tp, cfg.T_total);
        if isempty(intervals)
            template_scores(j) = 0;
            continue;
        end

        best_local = 0;
        best_iou = 0;
        best_ds = inf;
        best_dd = inf;

        for k = 1:size(intervals, 1)
            ts = intervals(k, 1);
            te = intervals(k, 2);
            td = max(1, te - ts + 1);

            iou = interval_iou(obs_start, obs_end, ts, te);
            ds = abs(obs_start - ts);
            dd = abs(obs_dur - td);

            s = cfg.alpha1 * iou + ...
                cfg.alpha2 * exp(-(ds.^2) / (2 * cfg.sigma_t.^2 + cfg.eps_val)) + ...
                cfg.alpha3 * exp(-(dd.^2) / (2 * cfg.sigma_d.^2 + cfg.eps_val));

            if s > best_local
                best_local = s;
                best_iou = iou;
                best_ds = ds;
                best_dd = dd;
            end
        end

        % 周期模板 phase 项
        s_phase = 1;
        if strcmpi(get_field_default(tp, 'mode', 'unknown'), 'periodic')
            dp = phase_delta(obs_start, tp);
            s_phase = exp(-(dp.^2) / (2 * cfg.sigma_p.^2 + cfg.eps_val));
        end

        score_j = best_local * s_phase;
        template_scores(j) = min(max(score_j, 0), 1);

        diag_info.iou_best(j) = best_iou;
        diag_info.delta_start_best(j) = best_ds;
        diag_info.delta_duration_best(j) = best_dd;
        diag_info.phase_score(j) = s_phase;
    end

    [best_score, idx_best] = max(template_scores);
    if ~isempty(idx_best) && idx_best >= 1 && idx_best <= N
        best_template_key = get_field_default(TemplateList(idx_best), 'template_key', '');
    end
    if N >= 2
        ss = sort(template_scores, 'descend');
        second_best_score = ss(2);
    end
end


%% ==================== 局部函数 ====================

function cfg = fill_time_match_defaults(Config)
    cfg.T_total = 500;
    cfg.alpha1 = 0.50;
    cfg.alpha2 = 0.25;
    cfg.alpha3 = 0.25;
    cfg.sigma_t = 6;
    cfg.sigma_d = 6;
    cfg.sigma_p = 4;
    cfg.eps_val = 1e-12;

    if isfield(Config, 'm0') && isfield(Config.m0, 'T_total')
        cfg.T_total = Config.m0.T_total;
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'time_match')
        tm = Config.m3.time_match;
        cfg.alpha1 = get_field_default(tm, 'alpha1', cfg.alpha1);
        cfg.alpha2 = get_field_default(tm, 'alpha2', cfg.alpha2);
        cfg.alpha3 = get_field_default(tm, 'alpha3', cfg.alpha3);
        cfg.sigma_t = get_field_default(tm, 'sigma_t', cfg.sigma_t);
        cfg.sigma_d = get_field_default(tm, 'sigma_d', cfg.sigma_d);
        cfg.sigma_p = get_field_default(tm, 'sigma_p', cfg.sigma_p);
    end
end


function [intervals, phase_delta_fn] = build_template_intervals(tp, T_total)
    intervals = zeros(0, 2);
    phase_delta_fn = @phase_delta_impl;

    mode = lower(get_field_default(tp, 'mode', 'unknown'));
    switch mode
        case 'periodic'
            period = scalarize_local(get_field_default(tp, 'period_frames', []));
            dur = scalarize_local(get_field_default(tp, 'duration_frames', []));
            phase = scalarize_local(get_field_default(tp, 'phase_frames', 0));
            if isempty(period) || isempty(dur) || period <= 0 || dur <= 0
                return;
            end
            starts = (1 + phase):period:T_total;
            for i = 1:numel(starts)
                s = starts(i);
                e = min(T_total, s + dur - 1);
                if e >= s
                    intervals(end+1, :) = [s, e]; %#ok<AGROW>
                end
            end

        case 'table'
            sch = get_field_default(tp, 'schedule', []);
            sch = unique(round(sch(:).'));
            sch = sch(sch >= 1 & sch <= T_total);
            if isempty(sch), return; end
            active = false(1, T_total);
            active(sch) = true;
            intervals = logical_to_segments(active);

        otherwise
            intervals = zeros(0, 2);
    end
end


function d = phase_delta_impl(obs_start, tp)
    period = scalarize_local(get_field_default(tp, 'period_frames', []));
    phase = scalarize_local(get_field_default(tp, 'phase_frames', 0));
    dur = scalarize_local(get_field_default(tp, 'duration_frames', 1));
    if isempty(period) || period <= 0
        d = 0;
        return;
    end
    p = mod(obs_start - 1 - phase, period);
    if p < dur
        d = 0;
    else
        d = min(abs(p - dur), abs(p - period));
    end
end


function iou = interval_iou(s1, e1, s2, e2)
    os = max(s1, s2);
    oe = min(e1, e2);
    if oe < os
        iou = 0;
        return;
    end
    inter = oe - os + 1;
    uni = max(e1, e2) - min(s1, s2) + 1;
    iou = inter / max(uni, 1);
end


function seg = logical_to_segments(active)
    seg = zeros(0, 2);
    T = numel(active);
    in_seg = false;
    s = 0;
    for t = 1:T
        if active(t) && ~in_seg
            in_seg = true;
            s = t;
        elseif ~active(t) && in_seg
            in_seg = false;
            seg(end+1, :) = [s, t - 1]; %#ok<AGROW>
        end
    end
    if in_seg
        seg(end+1, :) = [s, T]; %#ok<AGROW>
    end
end


function tp = resolve_time_prior_local(src)
    tp = struct();
    tp.mode = get_field_default(src, 'schedule_mode', 'unknown');
    tp.schedule = get_field_default(src, 'schedule_set', []);
    tp.period_frames = get_field_default(src, 'period_frames', []);
    tp.duration_frames = get_field_default(src, 'duration_frames', []);
    tp.phase_frames = get_field_default(src, 'phase_frames', []);

    if isfield(src, 'time_prior') && isstruct(src.time_prior)
        tp2 = src.time_prior;
        tp.mode = get_field_default(tp2, 'mode', tp.mode);
        tp.schedule = get_field_default(tp2, 'schedule', tp.schedule);
        tp.period_frames = get_field_default(tp2, 'period_frames', tp.period_frames);
        tp.duration_frames = get_field_default(tp2, 'duration_frames', tp.duration_frames);
        tp.phase_frames = get_field_default(tp2, 'phase_frames', tp.phase_frames);
    end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end


function v = scalarize_local(x)
    if isempty(x)
        v = [];
        return;
    end
    if isscalar(x)
        v = x;
    else
        v = x(1);
    end
end
