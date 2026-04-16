function [EventList, GroupList] = run_m3_event_classifier_single_source( ...
    Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config)
% RUN_M3_EVENT_CLASSIFIER_SINGLE_SOURCE  M3 入口（legacy + group）
%
%   新模式（grouping.enable=true）流程：
%     detect_band_events_m3
%       -> group_events_multiband_m3
%       -> extract_group_features_m3
%       -> compute_group_position_match_m3
%       -> compute_group_time_match_m3
%       -> classify_group_m3
%       -> route_group_m3

    fprintf('\n============================================\n');
    fprintf('  M3 事件分类模块\n');
    fprintf('============================================\n\n');

    tic;
    cfg = fill_run_defaults_m3(Config);

    % 1) 单 band 事件检测
    EventListRaw = detect_band_events_m3(Y_dBm_all, Y_lin_all, Config);
    if isempty(EventListRaw)
        fprintf('[M3] 未检测到事件\n');
        EventList = [];
        GroupList = [];
        return;
    end

    % 2) 单事件基础特征（底层）
    EventFeatList = extract_event_features_m3( ...
        EventListRaw, Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config);

    % legacy path
    if ~cfg.grouping_enable
        EventList = score_event_type_m3(EventFeatList, Config);
        EventList = route_event_m3(EventList, Config);
        GroupList = [];

        fprintf('\n[M3] legacy 模式完成: %d events, %.2f s\n', numel(EventList), toc);
        return;
    end

    % 3) group 级主流程
    GroupList = group_events_multiband_m3(EventFeatList, Config);
    GroupList = extract_group_features_m3( ...
        GroupList, EventFeatList, Y_dBm_all, Y_lin_all, SpatialFP, Config);
    GroupList = apply_group_template_matching_m3( ...
        GroupList, EventFeatList, SpatialFP, SignatureLib, SourceTemplates, Config);
    GroupList = classify_group_m3(GroupList, Config);
    GroupList = route_group_m3(GroupList, Config);

    % 4) 回填到单事件，保持 M4 兼容
    EventList = backfill_group_results_m3(EventFeatList, GroupList);

    elapsed = toc;
    fprintf('\n[M3] group 模式完成: %d groups, %d events, %.2f s\n', ...
        numel(GroupList), numel(EventList), elapsed);
end


%% ==================== 局部函数 ====================

function cfg = fill_run_defaults_m3(Config)
    cfg.grouping_enable = true;
    if isfield(Config, 'm3') && isfield(Config.m3, 'grouping') && ...
       isfield(Config.m3.grouping, 'enable')
        cfg.grouping_enable = logical(Config.m3.grouping.enable);
    end
end


function GroupList = apply_group_template_matching_m3( ...
    GroupList, EventList, SpatialFP, SignatureLib, SourceTemplates, Config)
% APPLY_GROUP_TEMPLATE_MATCHING_M3  逐 group 做位置/时间逐模板匹配

    if isempty(GroupList)
        return;
    end

    GroupList = ensure_group_fields_template_matching(GroupList);
    tpls = collect_template_sets_m3(SourceTemplates, SignatureLib);

    for g = 1:numel(GroupList)
        grp = GroupList(g);

        [sc_tr, key_tr, best_tr, second_tr, pos_diag_tr] = ...
            compute_group_position_match_m3(grp, EventList, SpatialFP, tpls.trusted, Config);
        [sc_pp, key_pp, best_pp, second_pp, pos_diag_pp] = ...
            compute_group_position_match_m3(grp, EventList, SpatialFP, tpls.prior_pos, Config);
        [~, ~, ~, ~, pos_diag_base] = ...
            compute_group_position_match_m3(grp, EventList, SpatialFP, struct([]), Config);

        [sc_tm, key_tm, best_tm, second_tm, time_diag] = ...
            compute_group_time_match_m3(grp, tpls.prior_time, Config);

        grp.trusted_template_keys = template_keys_local(tpls.trusted);
        grp.prior_pos_template_keys = template_keys_local(tpls.prior_pos);
        grp.prior_time_template_keys = template_keys_local(tpls.prior_time);

        grp.trusted_template_scores = sc_tr;
        grp.prior_pos_template_scores = sc_pp;
        grp.prior_time_template_scores = sc_tm;

        grp.trusted_best_template_key = key_tr;
        grp.prior_pos_best_template_key = key_pp;
        grp.prior_time_best_template_key = key_tm;
        grp.trusted_best_score = best_tr;
        grp.prior_pos_best_score = best_pp;
        grp.prior_time_best_score = best_tm;

        grp.time_template_scores = sc_tm;
        grp.position_template_scores = [sc_tr, sc_pp];

        if best_tr >= best_pp
            grp.best_position_template_key = key_tr;
            grp.best_position_score = best_tr;
            grp.second_best_position_score = max([second_tr, best_pp, 0]);
        else
            grp.best_position_template_key = key_pp;
            grp.best_position_score = best_pp;
            grp.second_best_position_score = max([second_pp, best_tr, 0]);
        end

        grp.best_time_template_key = key_tm;
        grp.best_time_score = best_tm;
        grp.second_best_time_score = second_tm;

        spread = NaN;
        if isfield(pos_diag_base, 'multiband_position_consistency')
            spread = pos_diag_base.multiband_position_consistency;
        end
        grp.multiband_position_spread = spread;
        grp.multiband_position_consistency = spread;

        grp.position_match_diag = pos_diag_base;
        grp.position_match_diag_trusted = pos_diag_tr;
        grp.position_match_diag_prior_pos = pos_diag_pp;
        grp.time_match_diag = time_diag;

        GroupList(g) = grp;
    end
end


function tpls = collect_template_sets_m3(SourceTemplates, SignatureLib)
% COLLECT_TEMPLATE_SETS_M3  收集三类模板（统一字段）

    tpls = struct();
    tpls.trusted = flatten_templates_from_source(SourceTemplates, 'trusted');
    tpls.prior_pos = flatten_templates_from_source(SourceTemplates, 'prior_pos');
    tpls.prior_time = flatten_templates_from_source(SourceTemplates, 'prior_time');

    if isstruct(SignatureLib) && isfield(SignatureLib, 'templates')
        st = SignatureLib.templates;
        if isempty(tpls.trusted)
            tpls.trusted = flatten_templates_from_signature(st, 'trusted_fixed', '');
        end
        if isempty(tpls.prior_pos)
            tpls.prior_pos = flatten_templates_from_signature(st, '', 'prior_pos_known');
        end
        if isempty(tpls.prior_time)
            tpls.prior_time = flatten_templates_from_signature(st, '', 'prior_time_known');
        end
    end
end


function out = flatten_templates_from_source(SourceTemplates, field_name)
    out = init_template_struct_array();
    if ~isfield(SourceTemplates, field_name) || isempty(SourceTemplates.(field_name))
        return;
    end
    arr = SourceTemplates.(field_name);
    for i = 1:numel(arr)
        out = append_template_local(out, arr(i));
    end
end


function out = flatten_templates_from_signature(st, type_name, subtype_name)
    out = init_template_struct_array();
    if isempty(st)
        return;
    end
    for i = 1:numel(st)
        keep = true;
        if ~isempty(type_name)
            keep = keep && strcmp(get_field_default(st(i), 'source_type', ''), type_name);
        end
        if ~isempty(subtype_name)
            keep = keep && strcmp(get_field_default(st(i), 'source_subtype', ''), subtype_name);
        end
        if keep
            out = append_template_local(out, st(i));
        end
    end
end


function out = init_template_struct_array()
    out = struct( ...
        'template_key', {}, ...
        'band_mask', {}, ...
        'position_prior', {}, ...
        'candidate_positions', {}, ...
        'time_prior', {}, ...
        'source_type', {}, ...
        'source_subtype', {});
end


function out = append_template_local(out, src)
    key = get_field_default(src, 'template_key', '');
    if isempty(key)
        return;
    end

    n = normalize_template_entry_local(src);
    k = numel(out) + 1;
    out(k).template_key = n.template_key;
    out(k).band_mask = n.band_mask;
    out(k).position_prior = n.position_prior;
    out(k).candidate_positions = n.candidate_positions;
    out(k).time_prior = n.time_prior;
    out(k).source_type = n.source_type;
    out(k).source_subtype = n.source_subtype;
end


function n = normalize_template_entry_local(src)
    n = struct();
    n.template_key = get_field_default(src, 'template_key', '');
    n.band_mask = get_band_mask_local(src);
    n.position_prior = get_position_prior_local(src);
    n.candidate_positions = n.position_prior.candidate_positions;
    n.time_prior = get_time_prior_local(src);
    n.source_type = get_field_default(src, 'source_type', '');
    n.source_subtype = get_field_default(src, 'source_subtype', '');
end


function bm = get_band_mask_local(src)
    bm = get_field_default(src, 'band_mask', []);
    if ~isempty(bm)
        bm = bm(:).';
        return;
    end
    b = get_field_default(src, 'band_id', []);
    if isempty(b)
        bm = [];
        return;
    end
    b = max(1, round(b(1)));
    bm = zeros(1, b);
    bm(b) = 1;
end


function pp = get_position_prior_local(src)
    pp = struct();
    pp.mode = 'unknown';
    pp.candidate_positions = [];
    pp.match_radius = [];

    if isfield(src, 'position_prior') && isstruct(src.position_prior)
        p = src.position_prior;
        pp.mode = get_field_default(p, 'mode', pp.mode);
        pp.candidate_positions = get_field_default(p, 'candidate_positions', pp.candidate_positions);
        pp.match_radius = get_field_default(p, 'match_radius', pp.match_radius);
        return;
    end

    pp.mode = get_field_default(src, 'known_position_mode', pp.mode);
    pp.candidate_positions = get_field_default(src, 'candidate_positions', pp.candidate_positions);
end


function tp = get_time_prior_local(src)
    tp = struct();
    tp.mode = 'unknown';
    tp.schedule = [];
    tp.period_frames = [];
    tp.duration_frames = [];
    tp.phase_frames = [];

    if isfield(src, 'time_prior') && isstruct(src.time_prior)
        t = src.time_prior;
        tp.mode = get_field_default(t, 'mode', tp.mode);
        tp.schedule = get_field_default(t, 'schedule', tp.schedule);
        tp.period_frames = get_field_default(t, 'period_frames', tp.period_frames);
        tp.duration_frames = get_field_default(t, 'duration_frames', tp.duration_frames);
        tp.phase_frames = get_field_default(t, 'phase_frames', tp.phase_frames);
        return;
    end

    tp.mode = get_field_default(src, 'schedule_mode', tp.mode);
    if strcmp(tp.mode, 'unknown') && strcmp(get_field_default(src, 'time_pattern_type', ''), 'scheduled')
        tp.mode = 'periodic';
    end
    tp.schedule = get_field_default(src, 'schedule_set', get_field_default(src, 'schedule', tp.schedule));
    tp.period_frames = get_field_default(src, 'period_frames', tp.period_frames);
    tp.duration_frames = get_field_default(src, 'duration_frames', tp.duration_frames);
    tp.phase_frames = get_field_default(src, 'phase_frames', tp.phase_frames);
end


function keys = template_keys_local(tpls)
    keys = {};
    if isempty(tpls)
        return;
    end
    keys = cell(1, numel(tpls));
    for i = 1:numel(tpls)
        keys{i} = get_field_default(tpls(i), 'template_key', '');
    end
end


function EventList = backfill_group_results_m3(EventFeatList, GroupList)
% BACKFILL_GROUP_RESULTS_M3  group 结果回填到事件列表

    EventList = EventFeatList;
    n_events = numel(EventList);

    for e = 1:n_events
        EventList(e).parent_group_id = 0;
        EventList(e).type_hat_group = 'unknown';
        EventList(e).linked_template_key_group = '';
        EventList(e).route_action_group = 'hold';
        EventList(e).hold_reason_group = 'no_group_assigned';

        EventList(e).type_hat = 'unknown';
        EventList(e).route_action = 'hold';
        EventList(e).linked_template_key = '';
        EventList(e).upgrade_hint = 'none';
        EventList(e).score_trusted = 0;
        EventList(e).score_prior_pos = 0;
        EventList(e).score_prior_time = 0;
        EventList(e).score_target = 0;

        EventList(e).group_id = 0;
        EventList(e).band_mask_group = [];
        EventList(e).active_band_count_group = 0;
        EventList(e).duration_frames_group = 0;
        EventList(e).power_excess_mean_dB_per_band = [];
        EventList(e).power_stability_per_band = [];
        EventList(e).n_valid_ap_per_band = [];
        EventList(e).power_zscore_per_band = [];
        EventList(e).multiband_position_spread = NaN;

        EventList(e).best_position_template_key = '';
        EventList(e).best_position_score = 0;
        EventList(e).second_best_position_score = 0;
        EventList(e).best_time_template_key = '';
        EventList(e).best_time_score = 0;
        EventList(e).second_best_time_score = 0;
        EventList(e).trusted_gate_flags = struct();
        EventList(e).trusted_subscores = struct();
        EventList(e).trusted_candidate_flag = false;
        EventList(e).trusted_hard_pass = false;
        EventList(e).trusted_hard_min_power_lin = NaN;
        EventList(e).trusted_hard_power_matrix_lin = [];
        EventList(e).trusted_hard_active_band_count = 0;
        EventList(e).trusted_hard_duration_frames = 0;
        EventList(e).trusted_hard_removed_from_scoring = false;
        EventList(e).pos_margin = 0;
        EventList(e).time_margin = 0;
        EventList(e).ratio_sum_lin_per_band = [];
        EventList(e).ratio_top1_lin_per_band = [];
        EventList(e).ratio_topK_lin_per_band = [];
        EventList(e).n_valid_ap_lin_per_band = [];
    end

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        idx_set = get_field_default(grp, 'member_event_indices', []);
        idx_set = idx_set(:).';
        idx_set = idx_set(idx_set >= 1 & idx_set <= n_events);
        for ii = 1:numel(idx_set)
            eidx = idx_set(ii);
            EventList(eidx).parent_group_id = grp.group_id;
            EventList(eidx).type_hat_group = grp.type_hat_group;
            EventList(eidx).linked_template_key_group = grp.linked_template_key;
            EventList(eidx).route_action_group = grp.route_action_group;
            EventList(eidx).hold_reason_group = grp.hold_reason;

            EventList(eidx).type_hat = grp.type_hat_group;
            EventList(eidx).route_action = grp.route_action_group;
            EventList(eidx).linked_template_key = grp.linked_template_key;
            EventList(eidx).score_trusted = get_field_default(grp, 'trusted_best_score', 0);
            EventList(eidx).score_prior_pos = get_field_default(grp, 'prior_pos_best_score', 0);
            EventList(eidx).score_prior_time = get_field_default(grp, 'prior_time_best_score', 0);
            EventList(eidx).score_target = double(strcmp(grp.type_hat_group, 'ordinary_target'));

            EventList(eidx).group_id = grp.group_id;
            EventList(eidx).band_mask_group = grp.band_mask;
            EventList(eidx).active_band_count_group = grp.active_band_count;
            EventList(eidx).duration_frames_group = grp.duration_frames;
            EventList(eidx).power_excess_mean_dB_per_band = grp.power_excess_mean_dB_per_band;
            EventList(eidx).power_stability_per_band = grp.power_stability_per_band;
            EventList(eidx).n_valid_ap_per_band = grp.n_valid_ap_per_band;
            EventList(eidx).power_zscore_per_band = get_field_default(grp, 'power_zscore_per_band', []);
            EventList(eidx).multiband_position_spread = get_field_default(grp, 'multiband_position_spread', NaN);

            EventList(eidx).best_position_template_key = get_field_default(grp, 'best_position_template_key', '');
            EventList(eidx).best_position_score = get_field_default(grp, 'best_position_score', 0);
            EventList(eidx).second_best_position_score = get_field_default(grp, 'second_best_position_score', 0);
            EventList(eidx).best_time_template_key = get_field_default(grp, 'best_time_template_key', '');
            EventList(eidx).best_time_score = get_field_default(grp, 'best_time_score', 0);
            EventList(eidx).second_best_time_score = get_field_default(grp, 'second_best_time_score', 0);
            EventList(eidx).trusted_gate_flags = get_field_default(grp, 'trusted_gate_flags', struct());
            EventList(eidx).trusted_subscores = get_field_default(grp, 'trusted_subscores', struct());
            EventList(eidx).trusted_candidate_flag = get_field_default(grp, 'trusted_candidate_flag', false);
            EventList(eidx).trusted_hard_pass = get_field_default(grp, 'trusted_hard_pass', false);
            EventList(eidx).trusted_hard_min_power_lin = get_field_default(grp, 'trusted_hard_min_power_lin', NaN);
            EventList(eidx).trusted_hard_power_matrix_lin = get_field_default(grp, 'trusted_hard_power_matrix_lin', []);
            EventList(eidx).trusted_hard_active_band_count = get_field_default(grp, 'trusted_hard_active_band_count', 0);
            EventList(eidx).trusted_hard_duration_frames = get_field_default(grp, 'trusted_hard_duration_frames', 0);
            EventList(eidx).trusted_hard_removed_from_scoring = ...
                get_field_default(grp, 'trusted_hard_removed_from_scoring', false);
            EventList(eidx).pos_margin = get_field_default(grp, 'pos_margin', 0);
            EventList(eidx).time_margin = get_field_default(grp, 'time_margin', 0);
            EventList(eidx).ratio_sum_lin_per_band = get_field_default(grp, 'ratio_sum_lin_per_band', []);
            EventList(eidx).ratio_top1_lin_per_band = get_field_default(grp, 'ratio_top1_lin_per_band', []);
            EventList(eidx).ratio_topK_lin_per_band = get_field_default(grp, 'ratio_topK_lin_per_band', []);
            EventList(eidx).n_valid_ap_lin_per_band = get_field_default(grp, 'n_valid_ap_lin_per_band', []);
        end
    end
end


function GroupList = ensure_group_fields_template_matching(GroupList)
    defs = struct();
    defs.trusted_template_keys = {};
    defs.prior_pos_template_keys = {};
    defs.prior_time_template_keys = {};
    defs.trusted_template_scores = [];
    defs.prior_pos_template_scores = [];
    defs.prior_time_template_scores = [];
    defs.position_template_scores = [];
    defs.time_template_scores = [];
    defs.trusted_best_template_key = '';
    defs.prior_pos_best_template_key = '';
    defs.prior_time_best_template_key = '';
    defs.trusted_best_score = 0;
    defs.prior_pos_best_score = 0;
    defs.prior_time_best_score = 0;
    defs.best_position_template_key = '';
    defs.best_position_score = 0;
    defs.second_best_position_score = 0;
    defs.best_time_template_key = '';
    defs.best_time_score = 0;
    defs.second_best_time_score = 0;
    defs.multiband_position_spread = NaN;
    defs.multiband_position_consistency = NaN;
    defs.position_match_diag = struct();
    defs.position_match_diag_trusted = struct();
    defs.position_match_diag_prior_pos = struct();
    defs.time_match_diag = struct();

    fn = fieldnames(defs);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(GroupList, f)
            [GroupList.(f)] = deal(defs.(f));
        end
    end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
