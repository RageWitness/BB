function SourceContext = bind_external_labels_to_events( ...
    LabelTable, FrameStates, Config)
% BIND_EXTERNAL_LABELS_TO_EVENTS  将外部标签绑定到帧/事件/信源
%
%   SourceContext = bind_external_labels_to_events(LabelTable, FrameStates, Config)
%
%   新框架下，源类型不再由 M3 内部判别。
%   本函数将外部输入的 LabelTable 与仿真帧推进逻辑联动。
%
%   输出：
%     SourceContext - struct，包含：
%       .label_table           - 原始标签表（传递）
%       .n_labels              - 标签总数
%       .per_frame(t)          - 每帧有效标签集合
%       .per_frame(t).labels   - 该帧有效的标签子集
%       .per_frame(t).n_active - 该帧有效标签数
%       .per_band_kind_map     - (T x B) cell，每格为 source_kind 或 ''
%       .per_band_uid_map      - (T x B) cell，每格为 source_uid 或 ''
%       .per_band_template_map - (T x B) cell，每格为 template_key 或 ''
%       .per_band_pos_map      - (T x B) cell，每格为 position_hint 或 []

    T = numel(FrameStates);
    B = 4;
    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        B = Config.m0.num_bands;
    end

    SourceContext = struct();
    SourceContext.label_table = LabelTable;
    SourceContext.n_labels = numel(LabelTable);

    % 初始化 per-frame
    per_frame = struct('labels', cell(1, T), 'n_active', num2cell(zeros(1, T)));

    % 初始化 per-band maps
    kind_map     = cell(T, B);
    uid_map      = cell(T, B);
    template_map = cell(T, B);
    pos_map      = cell(T, B);
    for t = 1:T
        for b = 1:B
            kind_map{t, b} = '';
            uid_map{t, b} = '';
            template_map{t, b} = '';
            pos_map{t, b} = [];
        end
    end

    if isempty(LabelTable)
        fprintf('[BindLabels] 无外部标签，SourceContext 全空\n');
        SourceContext.per_frame = per_frame;
        SourceContext.per_band_kind_map = kind_map;
        SourceContext.per_band_uid_map = uid_map;
        SourceContext.per_band_template_map = template_map;
        SourceContext.per_band_pos_map = pos_map;
        return;
    end

    % 逐标签展开到帧
    for i = 1:numel(LabelTable)
        lbl = LabelTable(i);
        sf = max(1, lbl.start_frame);
        ef = min(T, lbl.end_frame);
        b  = lbl.band_id;

        for t = sf:ef
            % 检查冲突：单源场景下每帧每 band 只有一个源
            if ~isempty(kind_map{t, b}) && ~isempty(kind_map{t, b})
                % 已有标签，保留优先级更高的
                old_kind = kind_map{t, b};
                new_kind = lbl.source_kind;
                if source_priority(new_kind) < source_priority(old_kind)
                    % 新标签优先级更高，替换
                    kind_map{t, b}     = new_kind;
                    uid_map{t, b}      = lbl.source_uid;
                    template_map{t, b} = lbl.template_key;
                    pos_map{t, b}      = lbl.position_hint;
                end
                % 否则保留旧的
            else
                kind_map{t, b}     = lbl.source_kind;
                uid_map{t, b}      = lbl.source_uid;
                template_map{t, b} = lbl.template_key;
                pos_map{t, b}      = lbl.position_hint;
            end
        end
    end

    % 汇总 per_frame
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
    SourceContext.per_band_kind_map = kind_map;
    SourceContext.per_band_uid_map = uid_map;
    SourceContext.per_band_template_map = template_map;
    SourceContext.per_band_pos_map = pos_map;

    % 统计
    n_covered = sum(~cellfun(@isempty, kind_map), 'all');
    fprintf('[BindLabels] 标签绑定完成: %d 条标签覆盖 %d 帧-band 格\n', ...
        numel(LabelTable), n_covered);
end


function p = source_priority(kind)
% SOURCE_PRIORITY  源优先级（数字越小优先级越高）
    switch kind
        case 'trusted',    p = 1;
        case 'prior_pos',  p = 2;
        case 'prior_time', p = 3;
        case 'ordinary',   p = 4;
        otherwise,         p = 99;
    end
end
