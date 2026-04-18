function LabelTable = ingest_external_source_labels(raw_labels, Config)
% INGEST_EXTERNAL_SOURCE_LABELS  外部源类型输入接口
%
%   LabelTable = ingest_external_source_labels(raw_labels, Config)
%
%   接收外部输入的源类型标签，校验后输出结构化标签表。
%   新框架下，源类型不再由 M3 内部判别，而是完全由外部输入。
%
%   输入：
%     raw_labels - struct array，每个元素至少包含：
%       .source_uid      - 外部源唯一标识（字符串或数字）
%       .band_id         - 频带编号
%       .source_kind     - 源类型: 'trusted' / 'prior_pos' / 'prior_time' / 'ordinary'
%       .start_frame     - 生效起始帧
%       .end_frame       - 生效结束帧
%     可选字段：
%       .position_hint   - [x, y] 若外部已知位置
%       .template_key    - 若外部已知模板键
%       .metadata        - 任意扩展字段
%
%     Config - 全局配置
%
%   输出：
%     LabelTable - 校验通过的结构化标签表（struct array）

    if nargin < 2, Config = struct(); end

    valid_kinds = {'trusted', 'prior_pos', 'prior_time', 'ordinary'};
    T_total = 500;
    B_total = 4;
    if isfield(Config, 'm0')
        if isfield(Config.m0, 'T_total'), T_total = Config.m0.T_total; end
        if isfield(Config.m0, 'num_bands'), B_total = Config.m0.num_bands; end
    end

    if isempty(raw_labels)
        fprintf('[ExternalLabels] 外部标签输入为空\n');
        LabelTable = init_empty_label_table();
        return;
    end

    n = numel(raw_labels);
    LabelTable = init_empty_label_table();
    n_valid = 0;
    n_reject = 0;

    for i = 1:n
        lbl = raw_labels(i);

        % --- 必填字段校验 ---
        [ok, reason] = validate_single_label(lbl, valid_kinds, B_total, T_total);
        if ~ok
            fprintf('[ExternalLabels] 标签 #%d 被拒: %s\n', i, reason);
            n_reject = n_reject + 1;
            continue;
        end

        % --- 构建标准化条目 ---
        entry = struct();
        entry.source_uid  = ensure_string(lbl.source_uid);
        entry.band_id     = lbl.band_id;
        entry.source_kind = lower(strtrim(lbl.source_kind));
        entry.start_frame = round(lbl.start_frame);
        entry.end_frame   = round(lbl.end_frame);

        if isfield(lbl, 'position_hint') && ~isempty(lbl.position_hint)
            entry.position_hint = lbl.position_hint(:)';
        else
            entry.position_hint = [];
        end

        if isfield(lbl, 'template_key') && ~isempty(lbl.template_key)
            entry.template_key = char(lbl.template_key);
        else
            entry.template_key = '';
        end

        if isfield(lbl, 'metadata')
            entry.metadata = lbl.metadata;
        else
            entry.metadata = struct();
        end

        n_valid = n_valid + 1;
        if n_valid == 1
            LabelTable = entry;
        else
            LabelTable(n_valid) = entry;
        end
    end

    fprintf('[ExternalLabels] 输入 %d 条, 接受 %d 条, 拒绝 %d 条\n', n, n_valid, n_reject);
end


%% ==================== 局部函数 ====================

function tbl = init_empty_label_table()
    tbl = struct( ...
        'source_uid',    {}, ...
        'band_id',       {}, ...
        'source_kind',   {}, ...
        'start_frame',   {}, ...
        'end_frame',     {}, ...
        'position_hint', {}, ...
        'template_key',  {}, ...
        'metadata',      {});
end


function [ok, reason] = validate_single_label(lbl, valid_kinds, B_total, T_total)
    ok = true;
    reason = '';

    if ~isfield(lbl, 'source_uid') || isempty(lbl.source_uid)
        ok = false; reason = 'missing source_uid'; return;
    end
    if ~isfield(lbl, 'band_id') || isempty(lbl.band_id)
        ok = false; reason = 'missing band_id'; return;
    end
    if ~isfield(lbl, 'source_kind') || isempty(lbl.source_kind)
        ok = false; reason = 'missing source_kind'; return;
    end
    if ~isfield(lbl, 'start_frame') || isempty(lbl.start_frame)
        ok = false; reason = 'missing start_frame'; return;
    end
    if ~isfield(lbl, 'end_frame') || isempty(lbl.end_frame)
        ok = false; reason = 'missing end_frame'; return;
    end

    kind = lower(strtrim(lbl.source_kind));
    if ~any(strcmp(kind, valid_kinds))
        ok = false; reason = sprintf('invalid source_kind: %s', kind); return;
    end

    bid = lbl.band_id;
    if ~isnumeric(bid) || bid < 1 || bid > B_total
        ok = false; reason = sprintf('band_id=%g out of range [1,%d]', bid, B_total); return;
    end

    sf = round(lbl.start_frame);
    ef = round(lbl.end_frame);
    if sf < 1 || ef > T_total || sf > ef
        ok = false; reason = sprintf('frame range [%d,%d] invalid (T=%d)', sf, ef, T_total); return;
    end
end


function s = ensure_string(v)
    if isnumeric(v)
        s = sprintf('src_%d', v);
    elseif ischar(v)
        s = v;
    elseif isstring(v)
        s = char(v);
    else
        s = 'unknown';
    end
end
