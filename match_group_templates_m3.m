function GroupList = match_group_templates_m3( ...
    GroupList, EventList, SpatialFP, SignatureLib, SourceTemplates, Config)
% MATCH_GROUP_TEMPLATES_M3  group 级逐模板匹配
%
%   输出：
%     trusted_template_scores / prior_pos_template_scores / prior_time_template_scores
%     trusted_best_template_key / prior_pos_best_template_key / prior_time_best_template_key
%     trusted_best_score / prior_pos_best_score / prior_time_best_score
%     best_position_score / best_time_score

    %#ok<INUSD>
    if isempty(GroupList)
        return;
    end

    GroupList = ensure_group_fields_match(GroupList);
    tpl = flatten_template_sets(SourceTemplates, SignatureLib);

    for g = 1:numel(GroupList)
        grp = GroupList(g);

        [trusted_scores, ~, ~, pos_diag_tr] = compute_group_position_match_m3( ...
            grp, EventList, SpatialFP, tpl.trusted, Config);
        [prior_pos_scores, ~, ~, pos_diag_pp] = compute_group_position_match_m3( ...
            grp, EventList, SpatialFP, tpl.prior_pos, Config);
        [~, ~, ~, pos_diag_base] = compute_group_position_match_m3( ...
            grp, EventList, SpatialFP, struct([]), Config);

        [trusted_best_score, trusted_best_key] = pick_best(trusted_scores, tpl.trusted);
        [prior_pos_best_score, prior_pos_best_key] = pick_best(prior_pos_scores, tpl.prior_pos);

        [prior_time_scores, prior_time_best_key, prior_time_best_score, time_diag] = ...
            compute_group_time_match_m3(grp, tpl.prior_time, Config);

        grp.trusted_template_keys = template_keys_cell(tpl.trusted);
        grp.prior_pos_template_keys = template_keys_cell(tpl.prior_pos);
        grp.prior_time_template_keys = template_keys_cell(tpl.prior_time);

        grp.trusted_template_scores = trusted_scores;
        grp.prior_pos_template_scores = prior_pos_scores;
        grp.prior_time_template_scores = prior_time_scores;

        grp.trusted_best_template_key = trusted_best_key;
        grp.prior_pos_best_template_key = prior_pos_best_key;
        grp.prior_time_best_template_key = prior_time_best_key;

        grp.trusted_best_score = trusted_best_score;
        grp.prior_pos_best_score = prior_pos_best_score;
        grp.prior_time_best_score = prior_time_best_score;

        grp.best_position_score = max([trusted_best_score, prior_pos_best_score, 0]);
        grp.best_time_score = prior_time_best_score;

        grp.position_match_score_per_template = [trusted_scores, prior_pos_scores];
        grp.time_pattern_match_score_per_template = prior_time_scores;

        grp.position_match_diag = pos_diag_base;
        grp.position_match_diag_trusted = pos_diag_tr;
        grp.position_match_diag_prior_pos = pos_diag_pp;
        grp.time_match_diag = time_diag;

        if isfield(pos_diag_base, 'multiband_position_consistency') && ...
           ~isnan(pos_diag_base.multiband_position_consistency)
            grp.multiband_position_consistency = pos_diag_base.multiband_position_consistency;
        end

        GroupList(g) = grp;
    end
end


%% ==================== 局部函数 ====================

function tpl = flatten_template_sets(SourceTemplates, SignatureLib)
    tpl = struct();
    tpl.trusted = flatten_trusted_templates(SourceTemplates);
    tpl.prior_pos = flatten_prior_pos_templates(SourceTemplates);
    tpl.prior_time = flatten_prior_time_templates(SourceTemplates);

    % 若 SourceTemplates 缺失，回退到 SignatureLib
    if isempty(tpl.trusted) || isempty(tpl.prior_pos) || isempty(tpl.prior_time)
        if isstruct(SignatureLib) && isfield(SignatureLib, 'templates')
            st = SignatureLib.templates;
            if isempty(tpl.trusted)
                tpl.trusted = flatten_signature_templates(st(strcmp({st.source_type}, 'trusted_fixed')));
            end
            if isempty(tpl.prior_pos)
                tpl.prior_pos = flatten_signature_templates(st(strcmp({st.source_subtype}, 'prior_pos_known')));
            end
            if isempty(tpl.prior_time)
                tpl.prior_time = flatten_signature_templates(st(strcmp({st.source_subtype}, 'prior_time_known')));
            end
        end
    end
end


function out = flatten_trusted_templates(SourceTemplates)
    out = init_template_array_local();
    if ~isfield(SourceTemplates, 'trusted') || isempty(SourceTemplates.trusted)
        return;
    end
    arr = SourceTemplates.trusted;
    for i = 1:numel(arr)
        out = append_template_entry_local(out, arr(i));
    end
end


function out = flatten_prior_pos_templates(SourceTemplates)
    out = init_template_array_local();
    if ~isfield(SourceTemplates, 'prior_pos') || isempty(SourceTemplates.prior_pos)
        return;
    end
    arr = SourceTemplates.prior_pos;
    for i = 1:numel(arr)
        out = append_template_entry_local(out, arr(i));
    end
end


function out = flatten_prior_time_templates(SourceTemplates)
    out = init_template_array_local();
    if ~isfield(SourceTemplates, 'prior_time') || isempty(SourceTemplates.prior_time)
        return;
    end
    arr = SourceTemplates.prior_time;
    for i = 1:numel(arr)
        out = append_template_entry_local(out, arr(i));
    end
end


function out = flatten_signature_templates(arr)
    out = init_template_array_local();
    for i = 1:numel(arr)
        out = append_template_entry_local(out, arr(i));
    end
end


function out = init_template_array_local()
% INIT_TEMPLATE_ARRAY_LOCAL  创建统一字段模板数组骨架
    out = struct( ...
        'template_key', {}, ...
        'band_mask', {}, ...
        'position_prior', {}, ...
        'candidate_positions', {}, ...
        'time_prior', {}, ...
        'source_type', {}, ...
        'source_subtype', {});
end


function out = append_template_entry_local(out, src)
% APPEND_TEMPLATE_ENTRY_LOCAL  逐字段填充，避免结构体异构赋值
    if ~isstruct(src) || ~isfield(src, 'template_key') || isempty(src.template_key)
        return;
    end

    n = normalize_template_entry(src);
    k = numel(out) + 1;
    out(k).template_key = n.template_key;
    out(k).band_mask = n.band_mask;
    out(k).position_prior = n.position_prior;
    out(k).candidate_positions = n.candidate_positions;
    out(k).time_prior = n.time_prior;
    out(k).source_type = n.source_type;
    out(k).source_subtype = n.source_subtype;
end


function tpl = normalize_template_entry(src)
% NORMALIZE_TEMPLATE_ENTRY  将任意模板条目投影为统一字段结构
    tpl = struct();
    tpl.template_key = get_field_default_local(src, 'template_key', '');
    tpl.band_mask = resolve_band_mask_local(src);
    tpl.position_prior = resolve_position_prior_local(src);
    tpl.candidate_positions = tpl.position_prior.candidate_positions;
    tpl.time_prior = resolve_time_prior_local(src);
    tpl.source_type = get_field_default_local(src, 'source_type', '');
    tpl.source_subtype = get_field_default_local(src, 'source_subtype', '');
end


function bm = resolve_band_mask_local(src)
    bm = [];
    if isfield(src, 'band_mask') && ~isempty(src.band_mask)
        bm = src.band_mask(:).';
        return;
    end
    if isfield(src, 'band_id') && ~isempty(src.band_id) && isnumeric(src.band_id)
        b = max(1, round(src.band_id(1)));
        bm = zeros(1, b);
        bm(b) = 1;
        return;
    end
end


function pp = resolve_position_prior_local(src)
    pp = struct();
    pp.mode = 'unknown';
    pp.candidate_positions = [];
    pp.match_radius = [];

    if isfield(src, 'position_prior') && isstruct(src.position_prior)
        p = src.position_prior;
        pp.mode = get_field_default_local(p, 'mode', pp.mode);
        pp.candidate_positions = get_field_default_local(p, 'candidate_positions', pp.candidate_positions);
        pp.match_radius = get_field_default_local(p, 'match_radius', pp.match_radius);
        return;
    end

    if isfield(src, 'known_position_mode')
        pp.mode = src.known_position_mode;
    end
    if isfield(src, 'candidate_positions')
        pp.candidate_positions = src.candidate_positions;
    end
end


function tp = resolve_time_prior_local(src)
    tp = struct();
    tp.mode = 'unknown';
    tp.schedule = [];
    tp.period_frames = [];
    tp.duration_frames = [];
    tp.phase_frames = [];

    if isfield(src, 'time_prior') && isstruct(src.time_prior)
        t = src.time_prior;
        tp.mode = get_field_default_local(t, 'mode', tp.mode);
        tp.schedule = get_field_default_local(t, 'schedule', tp.schedule);
        tp.period_frames = get_field_default_local(t, 'period_frames', tp.period_frames);
        tp.duration_frames = get_field_default_local(t, 'duration_frames', tp.duration_frames);
        tp.phase_frames = get_field_default_local(t, 'phase_frames', tp.phase_frames);
        return;
    end

    if isfield(src, 'schedule_mode')
        tp.mode = src.schedule_mode;
    elseif isfield(src, 'time_pattern_type') && strcmp(src.time_pattern_type, 'scheduled')
        tp.mode = 'periodic';
    end
    if isfield(src, 'schedule_set')
        tp.schedule = src.schedule_set;
    elseif isfield(src, 'schedule')
        tp.schedule = src.schedule;
    end
    if isfield(src, 'period_frames')
        tp.period_frames = src.period_frames;
    end
    if isfield(src, 'duration_frames')
        tp.duration_frames = src.duration_frames;
    end
    if isfield(src, 'phase_frames')
        tp.phase_frames = src.phase_frames;
    end
end


function v = get_field_default_local(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end


function [best_score, best_key] = pick_best(scores, templates)
    best_score = 0;
    best_key = '';
    if isempty(scores) || isempty(templates)
        return;
    end
    [best_score, idx] = max(scores);
    if idx >= 1 && idx <= numel(templates) && isfield(templates(idx), 'template_key')
        best_key = templates(idx).template_key;
    end
end


function keys = template_keys_cell(tpls)
    keys = {};
    if isempty(tpls)
        return;
    end
    keys = cell(1, numel(tpls));
    for i = 1:numel(tpls)
        if isfield(tpls(i), 'template_key')
            keys{i} = tpls(i).template_key;
        else
            keys{i} = '';
        end
    end
end


function GroupList = ensure_group_fields_match(GroupList)
% ENSURE_GROUP_FIELDS_MATCH  统一 group 字段，避免异构赋值
    defs = struct();
    defs.trusted_template_keys = {};
    defs.prior_pos_template_keys = {};
    defs.prior_time_template_keys = {};
    defs.trusted_template_scores = [];
    defs.prior_pos_template_scores = [];
    defs.prior_time_template_scores = [];
    defs.trusted_best_template_key = '';
    defs.prior_pos_best_template_key = '';
    defs.prior_time_best_template_key = '';
    defs.trusted_best_score = 0;
    defs.prior_pos_best_score = 0;
    defs.prior_time_best_score = 0;
    defs.best_position_score = 0;
    defs.best_time_score = 0;
    defs.position_match_score_per_template = [];
    defs.time_pattern_match_score_per_template = [];
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
