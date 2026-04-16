function GroupList = classify_group_m3(GroupList, Config)
% CLASSIFY_GROUP_M3  group 级顺序硬门限判决
%
%   顺序（严格串行，前级通过后不再判后级）：
%     1) trusted_fixed     — 三重一致：功率+稳定+空间
%     2) prior_pos_known   — 位置模板唯一性主导
%     3) prior_time_known  — 时间模板唯一性主导
%     4) ordinary_target   — 可检测存在 + 模板不支持
%     5) unknown / hold

    if isempty(GroupList)
        return;
    end

    GroupList = ensure_group_fields_classify(GroupList);
    cfg = fill_classify_defaults(Config);

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        active = find(grp.band_mask > 0);

        mean_power_z = nanmean_local(grp.power_zscore_per_band(active));
        mean_jitter  = nanmean_local(grp.power_stability_per_band(active));
        max_power_z  = nanmax_local(grp.power_zscore_per_band(active));
        max_valid_ap = nanmax_local(grp.n_valid_ap_per_band(active));
        spread       = get_field_default(grp, 'multiband_position_spread', NaN);

        pos_margin  = grp.best_position_score - grp.second_best_position_score;
        time_margin = grp.best_time_score     - grp.second_best_time_score;

        % ========== 1) trusted_fixed ==========
        % A. 功率一致
        f_active_bands = grp.active_band_count >= cfg.trusted.min_active_bands;
        f_power_z      = mean_power_z >= cfg.trusted.min_power_z;
        % B. 时间稳定
        f_jitter       = mean_jitter <= cfg.trusted.max_power_jitter_dB;
        % C. 空间固定
        f_pos          = grp.best_position_score >= cfg.trusted.min_position_score;
        f_spread       = ~isnan(spread) && (spread <= cfg.trusted.max_position_spread_m);

        trusted_pass = all([f_active_bands, f_power_z, f_jitter, f_pos, f_spread]);

        trusted_gate_flags = struct( ...
            'active_band_count',        f_active_bands, ...
            'mean_power_zscore',        f_power_z, ...
            'mean_power_jitter',        f_jitter, ...
            'best_position_score',      f_pos, ...
            'multiband_position_spread', f_spread);

        trusted_subscores = struct( ...
            'mean_power_z',   mean_power_z, ...
            'mean_jitter_dB', mean_jitter, ...
            'best_pos_score', grp.best_position_score, ...
            'spread_m',       spread, ...
            'active_bands',   grp.active_band_count);

        % ========== 2) prior_pos_known ==========
        prior_pos_pass = false;
        if ~trusted_pass
            prior_pos_pass = ...
                (grp.best_position_score >= cfg.prior_pos.min_position_score) && ...
                (pos_margin >= cfg.prior_pos.min_margin);
        end

        % ========== 3) prior_time_known ==========
        prior_time_pass = false;
        if ~trusted_pass && ~prior_pos_pass
            prior_time_pass = ...
                (grp.best_time_score >= cfg.prior_time.min_time_score) && ...
                (time_margin >= cfg.prior_time.min_margin);
        end

        % ========== 4) ordinary_target ==========
        ordinary_pass = false;
        if ~trusted_pass && ~prior_pos_pass && ~prior_time_pass
            % A. 可检测存在
            condA_detect = ...
                (max_power_z >= cfg.ordinary.min_detect_z) || ...
                (max_valid_ap >= cfg.ordinary.min_valid_ap);
            condA_duration = grp.duration_frames >= cfg.ordinary.min_duration_frames;

            % B. 模板不支持（前三级都没通过，且分数也不够门限）
            condB_no_template = ...
                (grp.best_position_score < cfg.prior_pos.min_position_score) && ...
                (grp.best_time_score     < cfg.prior_time.min_time_score);

            ordinary_pass = condA_detect && condA_duration && condB_no_template;
        end

        % ========== 最终类别 ==========
        if trusted_pass
            grp.type_hat_group = 'trusted_fixed';
            if ~isempty(grp.best_position_template_key)
                grp.linked_template_key = grp.best_position_template_key;
            else
                grp.linked_template_key = get_field_default(grp, 'trusted_best_template_key', '');
            end
            grp.hold_reason = '';
        elseif prior_pos_pass
            grp.type_hat_group = 'prior_pos_known';
            grp.linked_template_key = grp.best_position_template_key;
            grp.hold_reason = '';
        elseif prior_time_pass
            grp.type_hat_group = 'prior_time_known';
            grp.linked_template_key = grp.best_time_template_key;
            grp.hold_reason = '';
        elseif ordinary_pass
            grp.type_hat_group = 'ordinary_target';
            grp.linked_template_key = '';
            grp.hold_reason = '';
        else
            grp.type_hat_group = 'unknown';
            grp.linked_template_key = '';
            if ~exist('condA_detect', 'var')
                condA_detect = false; condA_duration = false; condB_no_template = false;
            end
            grp.hold_reason = build_hold_reason(grp, condA_detect, condA_duration, condB_no_template);
        end

        % ========== 诊断字段 ==========
        grp.trusted_gate_flags       = trusted_gate_flags;
        grp.trusted_subscores        = trusted_subscores;
        grp.trusted_gate_pass        = trusted_pass;
        grp.prior_pos_gate_pass      = prior_pos_pass;
        grp.prior_time_gate_pass     = prior_time_pass;
        grp.ordinary_gate_pass       = ordinary_pass;
        grp.pos_margin               = pos_margin;
        grp.time_margin              = time_margin;
        grp.mean_power_z_group       = mean_power_z;
        grp.mean_power_jitter_dB_group = mean_jitter;
        grp.max_power_z_group        = max_power_z;
        grp.max_valid_ap_group       = max_valid_ap;

        GroupList(g) = grp;
    end
end


%% ==================== 局部函数 ====================

function cfg = fill_classify_defaults(Config)
    cfg.trusted.min_active_bands     = 3;
    cfg.trusted.min_power_z          = 1.5;
    cfg.trusted.max_power_jitter_dB  = 3;
    cfg.trusted.min_position_score   = 0.70;
    cfg.trusted.max_position_spread_m = 10;

    cfg.prior_pos.min_position_score = 0.60;
    cfg.prior_pos.min_margin         = 0.10;

    cfg.prior_time.min_time_score    = 0.60;
    cfg.prior_time.min_margin        = 0.10;

    cfg.ordinary.min_detect_z        = 0.8;
    cfg.ordinary.min_duration_frames = 3;
    cfg.ordinary.min_valid_ap        = 2;

    if isfield(Config, 'm3')
        m3 = Config.m3;
        if isfield(m3, 'trusted')
            t = m3.trusted;
            cfg.trusted.min_active_bands = get_field_default(t, 'min_active_bands', cfg.trusted.min_active_bands);
            cfg.trusted.min_power_z = get_field_default(t, 'min_power_z', cfg.trusted.min_power_z);
            cfg.trusted.max_power_jitter_dB = get_field_default(t, 'max_power_jitter_dB', cfg.trusted.max_power_jitter_dB);
            cfg.trusted.min_position_score = get_field_default(t, 'min_position_score', cfg.trusted.min_position_score);
            cfg.trusted.max_position_spread_m = get_field_default(t, 'max_position_spread_m', ...
                get_field_default(t, 'max_position_spread', cfg.trusted.max_position_spread_m));
        end
        if isfield(m3, 'prior_pos')
            p = m3.prior_pos;
            cfg.prior_pos.min_position_score = get_field_default(p, 'min_position_score', cfg.prior_pos.min_position_score);
            cfg.prior_pos.min_margin = get_field_default(p, 'min_margin', cfg.prior_pos.min_margin);
        end
        if isfield(m3, 'prior_time')
            p = m3.prior_time;
            cfg.prior_time.min_time_score = get_field_default(p, 'min_time_score', cfg.prior_time.min_time_score);
            cfg.prior_time.min_margin = get_field_default(p, 'min_margin', cfg.prior_time.min_margin);
        end
        if isfield(m3, 'ordinary')
            o = m3.ordinary;
            cfg.ordinary.min_detect_z = get_field_default(o, 'min_detect_z', cfg.ordinary.min_detect_z);
            cfg.ordinary.min_duration_frames = get_field_default(o, 'min_duration_frames', cfg.ordinary.min_duration_frames);
            cfg.ordinary.min_valid_ap = get_field_default(o, 'min_valid_ap', cfg.ordinary.min_valid_ap);
        end
    end
end


function s = build_hold_reason(~, condA_detect, condA_duration, condB_no_template)
    if ~condA_detect
        s = 'signal_too_weak';
        return;
    end
    if ~condA_duration
        s = 'duration_too_short';
        return;
    end
    if ~condB_no_template
        s = 'template_partial_match';
        return;
    end
    s = 'unknown';
end


function GroupList = ensure_group_fields_classify(GroupList)
    defs = struct();
    defs.type_hat_group = 'unknown';
    defs.linked_template_key = '';
    defs.hold_reason = '';
    defs.trusted_gate_flags = struct();
    defs.trusted_subscores = struct();
    defs.trusted_gate_pass = false;
    defs.prior_pos_gate_pass = false;
    defs.prior_time_gate_pass = false;
    defs.ordinary_gate_pass = false;
    defs.pos_margin = 0;
    defs.time_margin = 0;
    defs.mean_power_z_group = NaN;
    defs.mean_power_jitter_dB_group = NaN;
    defs.max_power_z_group = NaN;
    defs.max_valid_ap_group = NaN;

    fn = fieldnames(defs);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(GroupList, f)
            [GroupList.(f)] = deal(defs.(f));
        end
    end
end


function y = nanmean_local(x)
    x = x(~isnan(x) & isfinite(x));
    if isempty(x), y = NaN; else, y = mean(x); end
end


function y = nanmax_local(x)
    x = x(~isnan(x) & isfinite(x));
    if isempty(x), y = -inf; else, y = max(x); end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
