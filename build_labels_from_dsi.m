function ExternalLabelsRaw = build_labels_from_dsi(DSI_all, Config)
% BUILD_LABELS_FROM_DSI  从 DrivenSourceInput 构建外部标签（主链路）
%
%   ExternalLabelsRaw = build_labels_from_dsi(DSI_all, Config)
%
%   以 DSI_all (T x B cell) 为唯一标签来源，保留完整 prior 信息。

    LC = source_label_constants();
    [T, B] = size(DSI_all);

    expose = isfield(Config, 'debug') && ...
             isfield(Config.debug, 'expose_true_source_state') && ...
             Config.debug.expose_true_source_state;

    ExternalLabelsRaw = struct( ...
        'source_uid',     {}, ...
        'band_id',        {}, ...
        'label',          {}, ...
        'start_frame',    {}, ...
        'end_frame',      {}, ...
        'position_hint',  {}, ...
        'template_key',   {}, ...
        'location_prior', {}, ...
        'power_prior',    {}, ...
        'metadata',       {});

    % --- 1. 收集所有有源 DSI 条目 ---
    entries = struct('t', {}, 'b', {}, 'label', {}, 'instance_id', {}, ...
                     'template_key', {}, 'dsi', {});
    for t = 1:T
        for b = 1:B
            dsi = DSI_all{t, b};
            if dsi.label == 0, continue; end

            e.t = t;
            e.b = b;
            e.label = dsi.label;

            if ~isempty(dsi.debug_info.instance_id) && ...
               isnumeric(dsi.debug_info.instance_id) && ...
               ~isnan(dsi.debug_info.instance_id)
                e.instance_id = dsi.debug_info.instance_id;
            else
                e.instance_id = NaN;
            end

            if ischar(dsi.debug_info.template_key) && ~isempty(dsi.debug_info.template_key)
                e.template_key = dsi.debug_info.template_key;
            else
                e.template_key = 'unknown_template';
            end

            e.dsi = dsi;

            if isempty(entries)
                entries = e;
            else
                entries(end+1) = e; %#ok<AGROW>
            end
        end
    end

    if isempty(entries)
        fprintf('[DSI->Labels] 无有源 DSI 条目\n');
        return;
    end

    % --- 2. 分组: label + band_id + instance_id + template_key ---
    n_entries = numel(entries);
    group_keys = cell(1, n_entries);
    for i = 1:n_entries
        e = entries(i);
        if ~isnan(e.instance_id)
            group_keys{i} = sprintf('%s__label%d__b%d__inst%d', ...
                e.template_key, e.label, e.b, e.instance_id);
        else
            group_keys{i} = sprintf('%s__label%d__b%d', ...
                e.template_key, e.label, e.b);
        end
    end

    [unique_groups, ~, ic] = unique(group_keys, 'stable');

    % --- 3. 逐组按连续帧切段 ---
    for g = 1:numel(unique_groups)
        idx = find(ic == g);
        frames = [entries(idx).t];
        [frames_sorted, ord] = sort(frames);
        idx = idx(ord);

        break_pos = [0, find(diff(frames_sorted) > 1), numel(frames_sorted)];

        for si = 1:numel(break_pos)-1
            seg_range = (break_pos(si)+1) : break_pos(si+1);
            seg_idx = idx(seg_range);

            ref = entries(seg_idx(1));
            dsi_ref = ref.dsi;

            lbl = struct();

            % source_uid
            if ~isnan(ref.instance_id)
                lbl.source_uid = sprintf('%s__label%d__b%d__inst%d__seg%d', ...
                    ref.template_key, ref.label, ref.b, ref.instance_id, si);
            else
                lbl.source_uid = sprintf('%s__label%d__b%d__seg%d', ...
                    ref.template_key, ref.label, ref.b, si);
            end

            lbl.band_id     = ref.b;
            lbl.label       = ref.label;
            lbl.start_frame = min(frames_sorted(seg_range));
            lbl.end_frame   = max(frames_sorted(seg_range));
            lbl.template_key = ref.template_key;

            % position_hint
            if strcmp(dsi_ref.location_prior.type, 'exact') && ...
               ~isempty(dsi_ref.location_prior.value)
                lbl.position_hint = dsi_ref.location_prior.value(:)';
            elseif ~isempty(dsi_ref.debug_info.true_position)
                lbl.position_hint = dsi_ref.debug_info.true_position(:)';
            else
                lbl.position_hint = [];
            end

            % prior 原样复制
            lbl.location_prior = dsi_ref.location_prior;
            lbl.power_prior    = dsi_ref.power_prior;

            % metadata
            md = struct();
            md.source_type_name = dsi_ref.meta.source_type_name;
            if LC.name_map.isKey(ref.label)
                md.label_name = LC.name_map(ref.label);
            else
                md.label_name = 'unknown';
            end
            if ~isnan(ref.instance_id)
                md.instance_id = ref.instance_id;
            else
                md.instance_id = 0;
            end
            md.template_key          = ref.template_key;
            md.band_id               = ref.b;
            md.band_list             = dsi_ref.band_list;
            md.n_frames              = numel(seg_idx);
            md.first_timestamp       = min(frames_sorted(seg_range));
            md.last_timestamp        = max(frames_sorted(seg_range));
            md.has_location_prior    = ~strcmp(dsi_ref.location_prior.type, 'none');
            md.has_power_prior       = ~strcmp(dsi_ref.power_prior.type, 'none');
            md.location_prior_type   = dsi_ref.location_prior.type;
            md.power_prior_type      = dsi_ref.power_prior.type;
            md.is_calibration_source = dsi_ref.meta.is_calibration_source;
            md.is_target_source      = dsi_ref.meta.is_target_source;
            md.is_opportunistic_source = dsi_ref.meta.is_opportunistic_source;

            if expose
                md.true_position    = dsi_ref.debug_info.true_position;
                md.true_tx_power    = dsi_ref.debug_info.true_tx_power;
                md.tx_power_by_band = dsi_ref.debug_info.tx_power_by_band;
            end

            lbl.metadata = md;

            if isempty(ExternalLabelsRaw)
                ExternalLabelsRaw = lbl;
            else
                ExternalLabelsRaw(end+1) = lbl; %#ok<AGROW>
            end
        end
    end

    fprintf('[DSI->Labels] 从 %d 个有源 DSI 条目生成 %d 条标签\n', ...
        n_entries, numel(ExternalLabelsRaw));

    % 统计
    labels_arr = [ExternalLabelsRaw.label];
    for lb = LC.ALL_LABELS
        cnt = sum(labels_arr == lb);
        if cnt > 0
            fprintf('  label=%d (%s): %d 条\n', lb, LC.name_map(lb), cnt);
        end
    end
end
