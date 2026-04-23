function SourceContext = bind_external_labels_to_events( ...
    LabelTable, FrameStates, Config)
% BIND_EXTERNAL_LABELS_TO_EVENTS  将外部标签绑定到帧（新标签体系）
%
%   新增 per_band_location_prior_map / per_band_power_prior_map /
%   per_band_metadata_map / per_band_label_name_map / per_band_instance_id_map

    LC = source_label_constants();

    T = numel(FrameStates);
    B = 4;
    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        B = Config.m0.num_bands;
    end

    SourceContext = struct();
    SourceContext.label_table = LabelTable;
    SourceContext.n_labels = numel(LabelTable);

    per_frame = struct('labels', cell(1, T), 'n_active', num2cell(zeros(1, T)));

    label_map    = zeros(T, B);
    uid_map      = cell(T, B);
    template_map = cell(T, B);
    pos_map      = cell(T, B);
    location_prior_map = cell(T, B);
    power_prior_map    = cell(T, B);
    metadata_map       = cell(T, B);
    label_name_map     = cell(T, B);
    instance_id_map    = zeros(T, B);

    default_lp = struct('type', 'none', 'value', []);
    default_pp = struct('type', 'none', 'value', []);

    for t = 1:T
        for b = 1:B
            uid_map{t, b}            = '';
            template_map{t, b}       = '';
            pos_map{t, b}            = [];
            location_prior_map{t, b} = default_lp;
            power_prior_map{t, b}    = default_pp;
            metadata_map{t, b}       = struct();
            label_name_map{t, b}     = '';
        end
    end

    if isempty(LabelTable)
        fprintf('[BindLabels] 无外部标签，SourceContext 全空\n');
        SourceContext.per_frame = per_frame;
        SourceContext.per_band_label_map    = label_map;
        SourceContext.per_band_uid_map      = uid_map;
        SourceContext.per_band_template_map = template_map;
        SourceContext.per_band_pos_map      = pos_map;
        SourceContext.per_band_location_prior_map = location_prior_map;
        SourceContext.per_band_power_prior_map    = power_prior_map;
        SourceContext.per_band_metadata_map       = metadata_map;
        SourceContext.per_band_label_name_map     = label_name_map;
        SourceContext.per_band_instance_id_map    = instance_id_map;
        return;
    end

    for i = 1:numel(LabelTable)
        lbl = LabelTable(i);
        sf = max(1, lbl.start_frame);
        ef = min(T, lbl.end_frame);
        b  = lbl.band_id;

        if LC.name_map.isKey(lbl.label)
            lname = LC.name_map(lbl.label);
        else
            lname = 'unknown';
        end

        if isfield(lbl, 'metadata') && isfield(lbl.metadata, 'instance_id')
            inst_id = lbl.metadata.instance_id;
        else
            inst_id = 0;
        end

        lp = default_lp;
        if isfield(lbl, 'location_prior') && isstruct(lbl.location_prior)
            lp = lbl.location_prior;
        end

        pp = default_pp;
        if isfield(lbl, 'power_prior') && isstruct(lbl.power_prior)
            pp = lbl.power_prior;
        end

        md = struct();
        if isfield(lbl, 'metadata')
            md = lbl.metadata;
        end

        for t = sf:ef
            if label_map(t, b) == 0
                label_map(t, b)            = lbl.label;
                uid_map{t, b}              = lbl.source_uid;
                template_map{t, b}         = lbl.template_key;
                pos_map{t, b}              = lbl.position_hint;
                location_prior_map{t, b}   = lp;
                power_prior_map{t, b}      = pp;
                metadata_map{t, b}         = md;
                label_name_map{t, b}       = lname;
                instance_id_map(t, b)      = inst_id;
            end
        end
    end

    for t = 1:T
        active_labels = [];
        for i = 1:numel(LabelTable)
            lbl = LabelTable(i);
            if t >= lbl.start_frame && t <= lbl.end_frame
                if isempty(active_labels)
                    active_labels = lbl;
                else
                    active_labels(end+1) = lbl; %#ok<AGROW>
                end
            end
        end
        if isempty(active_labels)
            per_frame(t).labels = struct([]);
        else
            per_frame(t).labels = active_labels;
        end
        per_frame(t).n_active = numel(active_labels);
    end

    SourceContext.per_frame = per_frame;
    SourceContext.per_band_label_map          = label_map;
    SourceContext.per_band_uid_map            = uid_map;
    SourceContext.per_band_template_map       = template_map;
    SourceContext.per_band_pos_map            = pos_map;
    SourceContext.per_band_location_prior_map = location_prior_map;
    SourceContext.per_band_power_prior_map    = power_prior_map;
    SourceContext.per_band_metadata_map       = metadata_map;
    SourceContext.per_band_label_name_map     = label_name_map;
    SourceContext.per_band_instance_id_map    = instance_id_map;

    n_covered = sum(label_map(:) > 0);
    fprintf('[BindLabels] 标签绑定完成: %d 条标签覆盖 %d 帧-band 格\n', ...
        numel(LabelTable), n_covered);
end
