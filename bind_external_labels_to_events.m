function SourceContext = bind_external_labels_to_events( ...
    LabelTable, FrameStates, Config)
% BIND_EXTERNAL_LABELS_TO_EVENTS  将外部标签绑定到帧（新标签体系）
%
%   SourceContext = bind_external_labels_to_events(LabelTable, FrameStates, Config)
%
%   使用数字标签 (1/2/3/4)，不再做优先级裁决（已在 M0 信源发生器中完成）。
%
%   输出：
%     SourceContext.label_table       - 原始标签表
%     SourceContext.n_labels          - 标签总数
%     SourceContext.per_frame(t).labels / .n_active
%     SourceContext.per_band_label_map    - (T x B) cell，数字标签或 0
%     SourceContext.per_band_uid_map      - (T x B) cell
%     SourceContext.per_band_template_map - (T x B) cell
%     SourceContext.per_band_pos_map      - (T x B) cell

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
    for t = 1:T
        for b = 1:B
            uid_map{t, b} = '';
            template_map{t, b} = '';
            pos_map{t, b} = [];
        end
    end

    if isempty(LabelTable)
        fprintf('[BindLabels] 无外部标签，SourceContext 全空\n');
        SourceContext.per_frame = per_frame;
        SourceContext.per_band_label_map = label_map;
        SourceContext.per_band_uid_map = uid_map;
        SourceContext.per_band_template_map = template_map;
        SourceContext.per_band_pos_map = pos_map;
        return;
    end

    for i = 1:numel(LabelTable)
        lbl = LabelTable(i);
        sf = max(1, lbl.start_frame);
        ef = min(T, lbl.end_frame);
        b  = lbl.band_id;

        for t = sf:ef
            if label_map(t, b) == 0
                label_map(t, b)    = lbl.label;
                uid_map{t, b}      = lbl.source_uid;
                template_map{t, b} = lbl.template_key;
                pos_map{t, b}      = lbl.position_hint;
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
    SourceContext.per_band_label_map = label_map;
    SourceContext.per_band_uid_map = uid_map;
    SourceContext.per_band_template_map = template_map;
    SourceContext.per_band_pos_map = pos_map;

    n_covered = sum(label_map(:) > 0);
    fprintf('[BindLabels] 标签绑定完成: %d 条标签覆盖 %d 帧-band 格\n', ...
        numel(LabelTable), n_covered);
end
