function GroupList = classify_group_m3(GroupList, Config)
% CLASSIFY_GROUP_M3  group-level ordered hard-gate classification

    if isempty(GroupList)
        return;
    end

    GroupList = ensure_group_fields_classify(GroupList);
    cfg = fill_classify_defaults(Config);

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        B = numel(grp.band_mask);
        active = find(grp.band_mask > 0);

        pos_margin  = grp.best_position_score - grp.second_best_position_score;
        time_margin = grp.best_time_score - grp.second_best_time_score;

        [trusted_hard_pass, trusted_diag] = eval_trusted_hard_gate(grp, cfg, B);
        grp.trusted_hard_pass = trusted_hard_pass;
        grp.trusted_hard_min_power_lin = trusted_diag.min_power_lin;
        grp.trusted_hard_active_band_count = trusted_diag.active_band_count;
        grp.trusted_hard_duration_frames = trusted_diag.duration_frames;

        if trusted_hard_pass
            grp.type_hat_group = 'trusted_fixed';
            grp.linked_template_key = choose_trusted_template_key(grp);
            grp.hold_reason = '';

            grp.trusted_hard_removed_from_scoring = true;
            grp.trusted_candidate_flag = true;  % backward-compatible alias
            grp.trusted_gate_pass = true;
            grp.prior_pos_gate_pass = false;
            grp.prior_time_gate_pass = false;
            grp.ordinary_gate_pass = false;

            grp.trusted_gate_flags = struct( ...
                'trusted_hard_pass', trusted_hard_pass, ...
                'all_bands_pass', trusted_diag.all_bands_pass, ...
                'all_aps_pass', trusted_diag.all_aps_pass, ...
                'duration_pass', trusted_diag.duration_pass, ...
                'mode', cfg.trusted_hard.mode, ...
                'aggregation', cfg.trusted_hard.aggregation);

            grp.trusted_subscores = struct( ...
                'trusted_hard_threshold_lin', cfg.trusted_hard.power_threshold_lin, ...
                'trusted_hard_min_power_lin', trusted_diag.min_power_lin, ...
                'trusted_hard_fraction_pass', trusted_diag.fraction_pass);

            grp.pos_margin = pos_margin;
            grp.time_margin = time_margin;
            GroupList(g) = grp;
            continue;
        end

        grp.trusted_hard_removed_from_scoring = false;
        grp.trusted_candidate_flag = false;

        prior_pos_pass = ...
            (grp.best_position_score >= cfg.prior_pos.min_position_score) && ...
            (pos_margin >= cfg.prior_pos.min_margin);

        prior_time_pass = false;
        if ~prior_pos_pass
            prior_time_pass = ...
                (grp.best_time_score >= cfg.prior_time.min_time_score) && ...
                (time_margin >= cfg.prior_time.min_margin);
        end

        ordinary_pass = false;
        condA_detect = false;
        condA_duration = false;
        condB_no_template = false;
        if ~prior_pos_pass && ~prior_time_pass
            r_sum = grp.ratio_sum_lin_per_band(active);
            r_top1 = grp.ratio_top1_lin_per_band(active);
            n_valid = grp.n_valid_ap_lin_per_band(active);

            condA_detect = ...
                any(r_sum  >= cfg.ordinary.min_ratio_sum) || ...
                any(r_top1 >= cfg.ordinary.min_ratio_top1) || ...
                any(n_valid >= cfg.ordinary.min_valid_ap);
            condA_duration = grp.duration_frames >= cfg.ordinary.min_duration_frames;

            condB_no_template = ...
                (grp.best_position_score < cfg.prior_pos.min_position_score) && ...
                (grp.best_time_score < cfg.prior_time.min_time_score);

            ordinary_pass = condA_detect && condA_duration && condB_no_template;
        end

        if prior_pos_pass
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
            grp.hold_reason = build_hold_reason(condA_detect, condA_duration, condB_no_template);
        end

        grp.trusted_gate_pass = false;
        grp.prior_pos_gate_pass = prior_pos_pass;
        grp.prior_time_gate_pass = prior_time_pass;
        grp.ordinary_gate_pass = ordinary_pass;

        grp.trusted_gate_flags = struct( ...
            'trusted_hard_pass', trusted_hard_pass, ...
            'all_bands_pass', trusted_diag.all_bands_pass, ...
            'all_aps_pass', trusted_diag.all_aps_pass, ...
            'duration_pass', trusted_diag.duration_pass, ...
            'mode', cfg.trusted_hard.mode, ...
            'aggregation', cfg.trusted_hard.aggregation);

        grp.trusted_subscores = struct( ...
            'trusted_hard_threshold_lin', cfg.trusted_hard.power_threshold_lin, ...
            'trusted_hard_min_power_lin', trusted_diag.min_power_lin, ...
            'trusted_hard_fraction_pass', trusted_diag.fraction_pass);

        grp.pos_margin = pos_margin;
        grp.time_margin = time_margin;
        GroupList(g) = grp;
    end
end


%% ==================== local helpers ====================

function [pass, diag] = eval_trusted_hard_gate(grp, cfg, B)
    diag = struct();
    diag.active_band_count = get_field_default(grp, 'active_band_count', 0);
    diag.duration_frames = get_field_default(grp, 'duration_frames', 0);
    diag.min_power_lin = get_field_default(grp, 'trusted_hard_min_power_lin', NaN);
    diag.fraction_pass = 0;
    diag.all_bands_pass = true;
    diag.all_aps_pass = false;
    diag.duration_pass = false;

    if ~cfg.trusted_hard.enable
        pass = false;
        return;
    end

    if cfg.trusted_hard.require_all_bands
        diag.all_bands_pass = (diag.active_band_count == B);
    else
        diag.all_bands_pass = (diag.active_band_count >= 1);
    end

    diag.duration_pass = diag.duration_frames >= cfg.trusted_hard.min_duration_frames;

    pwr = get_field_default(grp, 'trusted_hard_power_matrix_lin', []);
    if isempty(pwr)
        diag.all_aps_pass = false;
    else
        v = pwr(isfinite(pwr));
        if isempty(v)
            diag.all_aps_pass = false;
            diag.fraction_pass = 0;
            diag.min_power_lin = NaN;
        else
            diag.min_power_lin = min(v);
            pass_mask = v >= cfg.trusted_hard.power_threshold_lin;
            diag.fraction_pass = mean(double(pass_mask));

            switch lower(cfg.trusted_hard.mode)
                case 'fraction'
                    diag.all_aps_pass = (diag.fraction_pass >= cfg.trusted_hard.min_fraction);
                otherwise % strict
                    if cfg.trusted_hard.require_all_aps
                        diag.all_aps_pass = all(pass_mask);
                    else
                        diag.all_aps_pass = any(pass_mask);
                    end
            end
        end
    end

    pass = diag.all_bands_pass && diag.all_aps_pass && diag.duration_pass;
end


function key = choose_trusted_template_key(grp)
    key = '';
    key_pos = get_field_default(grp, 'best_position_template_key', '');
    key_tr = get_field_default(grp, 'trusted_best_template_key', '');
    if ~isempty(key_pos)
        key = key_pos;
    elseif ~isempty(key_tr)
        key = key_tr;
    end
end


function cfg = fill_classify_defaults(Config)
    cfg.trusted_hard.enable = true;
    cfg.trusted_hard.power_threshold_lin = 1e5;
    cfg.trusted_hard.require_all_bands = true;
    cfg.trusted_hard.require_all_aps = true;
    cfg.trusted_hard.min_duration_frames = 2;
    cfg.trusted_hard.aggregation = 'mean';
    cfg.trusted_hard.mode = 'strict';
    cfg.trusted_hard.min_fraction = 0.95;

    cfg.prior_pos.min_position_score = 0.60;
    cfg.prior_pos.min_margin = 0.10;

    cfg.prior_time.min_time_score = 0.60;
    cfg.prior_time.min_margin = 0.10;

    cfg.ordinary.min_ratio_sum = 3;
    cfg.ordinary.min_ratio_top1 = 5;
    cfg.ordinary.min_valid_ap = 2;
    cfg.ordinary.min_duration_frames = 3;

    if isfield(Config, 'm3')
        m3 = Config.m3;
        if isfield(m3, 'trusted_hard')
            t = m3.trusted_hard;
            cfg.trusted_hard.enable = gfd(t, 'enable', cfg.trusted_hard.enable);
            cfg.trusted_hard.power_threshold_lin = gfd(t, 'power_threshold_lin', cfg.trusted_hard.power_threshold_lin);
            cfg.trusted_hard.require_all_bands = gfd(t, 'require_all_bands', cfg.trusted_hard.require_all_bands);
            cfg.trusted_hard.require_all_aps = gfd(t, 'require_all_aps', cfg.trusted_hard.require_all_aps);
            cfg.trusted_hard.min_duration_frames = gfd(t, 'min_duration_frames', cfg.trusted_hard.min_duration_frames);
            cfg.trusted_hard.aggregation = gfd(t, 'aggregation', cfg.trusted_hard.aggregation);
            cfg.trusted_hard.mode = gfd(t, 'mode', cfg.trusted_hard.mode);
            cfg.trusted_hard.min_fraction = gfd(t, 'min_fraction', cfg.trusted_hard.min_fraction);
            cfg.trusted_hard.min_fraction = gfd(t, 'rho', cfg.trusted_hard.min_fraction);
        end

        if isfield(m3, 'prior_pos')
            p = m3.prior_pos;
            cfg.prior_pos.min_position_score = gfd(p, 'min_position_score', cfg.prior_pos.min_position_score);
            cfg.prior_pos.min_margin = gfd(p, 'min_margin', cfg.prior_pos.min_margin);
        end
        if isfield(m3, 'prior_time')
            p = m3.prior_time;
            cfg.prior_time.min_time_score = gfd(p, 'min_time_score', cfg.prior_time.min_time_score);
            cfg.prior_time.min_margin = gfd(p, 'min_margin', cfg.prior_time.min_margin);
        end
        if isfield(m3, 'ordinary')
            o = m3.ordinary;
            cfg.ordinary.min_ratio_sum = gfd(o, 'min_ratio_sum', cfg.ordinary.min_ratio_sum);
            cfg.ordinary.min_ratio_top1 = gfd(o, 'min_ratio_top1', cfg.ordinary.min_ratio_top1);
            cfg.ordinary.min_valid_ap = gfd(o, 'min_valid_ap', cfg.ordinary.min_valid_ap);
            cfg.ordinary.min_duration_frames = gfd(o, 'min_duration_frames', cfg.ordinary.min_duration_frames);
        end
    end
end


function s = build_hold_reason(condA_detect, condA_duration, condB_no_template)
    if ~condA_detect
        s = 'signal_too_weak';
    elseif ~condA_duration
        s = 'duration_too_short';
    elseif ~condB_no_template
        s = 'template_partial_match';
    else
        s = 'unknown';
    end
end


function GroupList = ensure_group_fields_classify(GroupList)
    defs = struct();
    defs.type_hat_group = 'unknown';
    defs.linked_template_key = '';
    defs.hold_reason = '';
    defs.trusted_candidate_flag = false;
    defs.trusted_gate_flags = struct();
    defs.trusted_subscores = struct();
    defs.trusted_gate_pass = false;
    defs.prior_pos_gate_pass = false;
    defs.prior_time_gate_pass = false;
    defs.ordinary_gate_pass = false;
    defs.pos_margin = 0;
    defs.time_margin = 0;

    defs.trusted_hard_pass = false;
    defs.trusted_hard_min_power_lin = NaN;
    defs.trusted_hard_power_matrix_lin = [];
    defs.trusted_hard_active_band_count = 0;
    defs.trusted_hard_duration_frames = 0;
    defs.trusted_hard_removed_from_scoring = false;

    fn = fieldnames(defs);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(GroupList, f)
            [GroupList.(f)] = deal(defs.(f));
        end
    end
end


function v = gfd(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
