function EventList = score_event_type_m3(EventList, Config)
% SCORE_EVENT_TYPE_M3  对每个事件计算四类源评分并判定类型
%
%   关键改动：
%     1) trusted_fixed 改为硬门限类（不再并行软竞争）
%     2) ordinary_target 改为兜底类（前三类失败后才可进入）
%     3) relative power 仅作为辅助证据，不直接定义类别

    sc = fill_score_defaults(Config);
    n_events = numel(EventList);

    % 预分配新字段
    for e = 1:n_events
        EventList(e).score_trusted = 0;
        EventList(e).score_prior_pos = 0;
        EventList(e).score_prior_time = 0;
        EventList(e).score_target = 0;
        EventList(e).type_hat = '';

        % 诊断字段
        EventList(e).active_band_count = 0;
        EventList(e).mean_power_excess_dB_ctx = NaN;
        EventList(e).mean_power_jitter_dB_ctx = NaN;
        EventList(e).position_spread_m_ctx = NaN;
        EventList(e).trusted_gate_pass = false;
        EventList(e).prior_pos_gate_pass = false;
        EventList(e).prior_time_gate_pass = false;
        EventList(e).ordinary_gate_pass = false;
    end

    for e = 1:n_events
        ev = EventList(e);

        % ===== 0) 上下文统计（跨重叠事件/频带）=====
        ctx = build_event_context(ev, EventList, sc);
        ev.active_band_count = ctx.active_band_count;
        ev.mean_power_excess_dB_ctx = ctx.mean_power_excess_dB;
        ev.mean_power_jitter_dB_ctx = ctx.mean_power_jitter_dB;
        ev.position_spread_m_ctx = ctx.position_spread_m;

        % ===== 1) 相对功率只做“信号显著性”辅助 =====
        signal_ok = is_signal_significant(ev, ctx, sc);
        stable_ok = ctx.mean_power_jitter_dB <= sc.trusted.max_power_jitter_dB;
        pos_ok = ev.position_prior_match_score >= sc.trusted.min_position_score;
        spread_ok = isnan(ctx.position_spread_m) || ...
                    (ctx.position_spread_m <= sc.trusted.max_position_spread_m);

        % ===== 2) trusted_fixed 硬门限 =====
        trusted_pass = ...
            (ctx.active_band_count >= sc.trusted.min_active_bands) && ...
            (ctx.mean_power_excess_dB >= sc.trusted.min_power_excess_dB) && ...
            stable_ok && ...
            pos_ok && ...
            spread_ok && ...
            signal_ok;

        % ===== 3) prior gating =====
        prior_pos_pass = ...
            (~trusted_pass) && ...
            signal_ok && ...
            (ev.schedule_match_score >= sc.prior_pos.min_schedule_score) && ...
            (ev.position_prior_match_score >= sc.prior_pos.min_position_score) && ...
            (ev.power_stability_est <= sc.prior_pos.max_power_jitter_dB);

        prior_time_pass = ...
            (~trusted_pass) && (~prior_pos_pass) && ...
            signal_ok && ...
            (ev.schedule_match_score >= sc.prior_time.min_schedule_score) && ...
            (ev.power_stability_est <= sc.prior_time.max_power_jitter_dB);

        % ===== 4) ordinary 兜底 gating =====
        ordinary_pass = ...
            (~trusted_pass) && (~prior_pos_pass) && (~prior_time_pass) && ...
            (ev.duration >= sc.ordinary.min_duration_frames) && ...
            (get_field_default(ev, 'n_valid_ap', 0) >= sc.ordinary.min_valid_ap) && ...
            signal_ok;

        % 记录 gate
        ev.trusted_gate_pass = trusted_pass;
        ev.prior_pos_gate_pass = prior_pos_pass;
        ev.prior_time_gate_pass = prior_time_pass;
        ev.ordinary_gate_pass = ordinary_pass;

        % ===== 5) 辅助评分（仅供 route/诊断，不做并行抢类）=====
        [score_tr, score_pp, score_pt, score_tg] = ...
            build_aux_scores(ev, ctx, sc, trusted_pass, prior_pos_pass, prior_time_pass, ordinary_pass, signal_ok);
        ev.score_trusted = score_tr;
        ev.score_prior_pos = score_pp;
        ev.score_prior_time = score_pt;
        ev.score_target = score_tg;

        % ===== 6) 分类顺序（硬判决）=====
        if trusted_pass
            ev.type_hat = 'trusted_fixed';
        elseif prior_pos_pass
            ev.type_hat = 'prior_pos_known';
        elseif prior_time_pass
            ev.type_hat = 'prior_time_known';
        elseif ordinary_pass
            ev.type_hat = 'ordinary_target';
        else
            % 不满足兜底条件时，不让 ordinary 误接收。
            % 保持 type 为 ordinary_target，但分数压低，route 阶段会进 hold。
            ev.type_hat = 'ordinary_target';
            ev.score_trusted = min(ev.score_trusted, sc.score_cap_fail);
            ev.score_prior_pos = min(ev.score_prior_pos, sc.score_cap_fail);
            ev.score_prior_time = min(ev.score_prior_time, sc.score_cap_fail);
            ev.score_target = min(ev.score_target, sc.score_cap_fail);
        end

        EventList(e) = ev;
    end

    if n_events > 0
        types = {EventList.type_hat};
        fprintf('[M3] 事件分类完成: trusted=%d, prior_pos=%d, prior_time=%d, target=%d\n', ...
            sum(strcmp(types, 'trusted_fixed')), ...
            sum(strcmp(types, 'prior_pos_known')), ...
            sum(strcmp(types, 'prior_time_known')), ...
            sum(strcmp(types, 'ordinary_target')));
    end
end


%% ==================== 局部函数 ====================

function ctx = build_event_context(ev, EventList, sc)
% BUILD_EVENT_CONTEXT  构建事件上下文统计（跨重叠事件）

    idx_ov = find_overlapped_events(ev, EventList);
    if isempty(idx_ov)
        idx_ov = 1; %#ok<NASGU>
        ov = ev;
    else
        ov = EventList(idx_ov);
    end

    % active band count：取 overlap 统计与 band_coverage 的较大值
    if ~isempty(ov)
        ctx.active_band_count = numel(unique([ov.band_id]));
    else
        ctx.active_band_count = sum(ev.band_coverage_vec > 0);
    end
    ctx.active_band_count = max(ctx.active_band_count, sum(ev.band_coverage_vec > 0));

    ex = arrayfun(@(x)get_field_default(x, 'power_excess_mean_dB', NaN), ov);
    ex = ex(~isnan(ex));
    if isempty(ex)
        ctx.mean_power_excess_dB = get_field_default(ev, 'power_excess_mean_dB', -inf);
    else
        ctx.mean_power_excess_dB = mean(ex);
    end

    jit = arrayfun(@(x)get_field_default(x, 'power_stability_est', NaN), ov);
    jit = jit(~isnan(jit));
    if isempty(jit)
        ctx.mean_power_jitter_dB = get_field_default(ev, 'power_stability_est', inf);
    else
        ctx.mean_power_jitter_dB = mean(jit);
    end

    % 若可获得多 band 位置结果，计算离散度；否则 NaN（自动放行）
    pos_xy = collect_positions_from_events(ov);
    if size(pos_xy, 1) >= 2
        c = mean(pos_xy, 1);
        d = sqrt(sum((pos_xy - c).^2, 2));
        ctx.position_spread_m = max(d);  % 最大离中心偏移
    else
        ctx.position_spread_m = NaN;
    end

    % 用于辅助评分
    if sc.total_bands > 0
        ctx.coverage_ratio = min(1, ctx.active_band_count / sc.total_bands);
    else
        ctx.coverage_ratio = 0;
    end
end


function idx = find_overlapped_events(ev, EventList)
% FIND_OVERLAPPED_EVENTS  查找时间上与 ev 重叠的事件
    n = numel(EventList);
    hit = false(1, n);
    for k = 1:n
        other = EventList(k);
        t0 = max(ev.t_start, other.t_start);
        t1 = min(ev.t_end, other.t_end);
        hit(k) = (t1 >= t0);
    end
    idx = find(hit);
end


function pos_xy = collect_positions_from_events(EventArray)
% COLLECT_POSITIONS_FROM_EVENTS  从事件中收集可用位置
%   目前 M3 阶段通常没有多 band 位置估计；若无则返回空。

    pos_xy = zeros(0, 2);
    for i = 1:numel(EventArray)
        e = EventArray(i);
        if isfield(e, 'est_pos_xy') && numel(e.est_pos_xy) == 2
            pos_xy(end+1, :) = reshape(e.est_pos_xy, 1, 2); %#ok<AGROW>
        elseif isfield(e, 'pos_est_xy') && numel(e.pos_est_xy) == 2
            pos_xy(end+1, :) = reshape(e.pos_est_xy, 1, 2); %#ok<AGROW>
        end
    end
end


function ok = is_signal_significant(ev, ctx, sc)
% IS_SIGNAL_SIGNIFICANT  relative power 作为“是否有有效信号”的辅助判据
    if sc.use_relative_power
        n_valid = get_field_default(ev, 'n_valid_ap', 0);
        ex_top1 = get_field_default(ev, 'power_excess_top1_dB', -inf);
        ex_mean = ctx.mean_power_excess_dB;
        ok = (n_valid >= sc.signal.min_valid_ap) && ...
             (ex_top1 >= sc.signal.min_top1_excess_dB) && ...
             (ex_mean >= sc.signal.min_mean_excess_dB);
    else
        ok = (ev.power_level_est >= sc.signal.min_abs_power_dBm);
    end
end


function [score_tr, score_pp, score_pt, score_tg] = build_aux_scores( ...
    ev, ctx, sc, trusted_pass, prior_pos_pass, prior_time_pass, ordinary_pass, signal_ok)
% BUILD_AUX_SCORES  构建辅助分数（不参与并行抢类）

    S_cover = ctx.coverage_ratio;
    S_time  = get_field_default(ev, 'schedule_match_score', 0);
    S_pos   = get_field_default(ev, 'position_prior_match_score', 0);
    S_stable = exp(-ctx.mean_power_jitter_dB / max(sc.trusted.max_power_jitter_dB, 1e-6));

    if sc.use_relative_power
        ex_mean = ctx.mean_power_excess_dB;
        ex_top1 = get_field_default(ev, 'power_excess_top1_dB', -inf);
        S_sig = 0.5 * sigmoid((ex_mean - sc.trusted.min_power_excess_dB) / sc.signal.excess_scale_dB) + ...
                0.5 * sigmoid((ex_top1 - sc.signal.min_top1_excess_dB) / sc.signal.excess_scale_dB);
    else
        S_sig = sigmoid((ev.power_level_est - sc.signal.min_abs_power_dBm) / 5);
    end

    S_dur = min(1, ev.duration / max(sc.ordinary.min_duration_frames, 1));

    tr_aux = clamp01(0.35*S_cover + 0.25*S_stable + 0.25*S_pos + 0.15*S_sig);
    pp_aux = clamp01(0.45*S_time + 0.35*S_pos + 0.20*S_stable);
    pt_aux = clamp01(0.60*S_time + 0.40*S_stable);
    tg_aux = clamp01(0.50*S_sig + 0.50*S_dur);

    if trusted_pass
        score_tr = 0.90 + 0.08 * tr_aux;
    else
        score_tr = 0.05 * tr_aux;
    end

    if prior_pos_pass
        score_pp = 0.82 + 0.10 * pp_aux;
    else
        score_pp = 0.05 * pp_aux;
    end

    if prior_time_pass
        score_pt = 0.78 + 0.10 * pt_aux;
    else
        score_pt = 0.05 * pt_aux;
    end

    if ordinary_pass
        score_tg = 0.72 + 0.12 * tg_aux;
    elseif signal_ok
        score_tg = 0.05 * tg_aux;
    else
        score_tg = 0.01 * tg_aux;
    end

    % 上限保护
    score_tr = min(score_tr, 0.99);
    score_pp = min(score_pp, 0.99);
    score_pt = min(score_pt, 0.99);
    score_tg = min(score_tg, 0.99);
end


function sc = fill_score_defaults(Config)
% FILL_SCORE_DEFAULTS  填充评分参数默认值

    if isfield(Config, 'm3') && isfield(Config.m3, 'score')
        s = Config.m3.score;
    else
        s = struct();
    end

    if isfield(Config, 'm3') && isfield(Config.m3, 'relative_power')
        rp = Config.m3.relative_power;
    else
        rp = struct();
    end
    if isfield(rp, 'enable')
        sc.use_relative_power = rp.enable;
    else
        sc.use_relative_power = true;
    end

    % band 数用于 coverage ratio
    sc.total_bands = 4;
    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        sc.total_bands = Config.m0.num_bands;
    end

    % ===== trusted 硬门限 =====
    if isfield(Config, 'm3') && isfield(Config.m3, 'trusted')
        t = Config.m3.trusted;
    else
        t = struct();
    end
    sc.trusted.min_active_bands      = get_field_default(t, 'min_active_bands', 3);
    sc.trusted.min_power_excess_dB   = get_field_default(t, 'min_power_excess_dB', 8);
    sc.trusted.max_power_jitter_dB   = get_field_default(t, 'max_power_jitter_dB', 3);
    sc.trusted.min_position_score    = get_field_default(t, 'min_position_score', 0.7);
    sc.trusted.max_position_spread_m = get_field_default(t, 'max_position_spread', 15);

    % ===== prior gate =====
    if isfield(Config, 'm3') && isfield(Config.m3, 'prior_pos')
        pp = Config.m3.prior_pos;
    else
        pp = struct();
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'prior_time')
        pt = Config.m3.prior_time;
    else
        pt = struct();
    end
    sc.prior_pos.min_schedule_score = get_field_default(pp, 'min_schedule_score', 0.45);
    sc.prior_pos.min_position_score = get_field_default(pp, 'min_position_score', 0.55);
    sc.prior_pos.max_power_jitter_dB = get_field_default(pp, 'max_power_jitter_dB', 6);

    sc.prior_time.min_schedule_score = get_field_default(pt, 'min_schedule_score', 0.55);
    sc.prior_time.max_power_jitter_dB = get_field_default(pt, 'max_power_jitter_dB', 6);

    % ===== ordinary 兜底 =====
    if isfield(Config, 'm3') && isfield(Config.m3, 'ordinary')
        og = Config.m3.ordinary;
    else
        og = struct();
    end
    sc.ordinary.min_duration_frames = get_field_default(og, 'min_duration_frames', 3);
    sc.ordinary.min_valid_ap = get_field_default(og, 'min_valid_ap', 2);

    % ===== signal 辅助判据 =====
    if isfield(Config, 'm3') && isfield(Config.m3, 'signal')
        sg = Config.m3.signal;
    else
        sg = struct();
    end
    sc.signal.min_valid_ap = get_field_default(sg, 'min_valid_ap', 1);
    sc.signal.min_top1_excess_dB = get_field_default(sg, 'min_top1_excess_dB', 3);
    sc.signal.min_mean_excess_dB = get_field_default(sg, 'min_mean_excess_dB', 1.5);
    sc.signal.min_abs_power_dBm = get_field_default(sg, 'min_abs_power_dBm', -110);
    sc.signal.excess_scale_dB = get_field_default(sg, 'excess_scale_dB', 3);

    % route hold 保护：分类失败时压低分数到此值以下
    sc.score_cap_fail = 0.05;

    %#ok<NASGU>
    s = s; % 保留，避免未来扩展误删
end


function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end


function y = clamp01(x)
    y = min(max(x, 0), 1);
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
