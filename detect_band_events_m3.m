function EventListRaw = detect_band_events_m3(Y_dBm_all, Y_lin_all, Config)
% DETECT_BAND_EVENTS_M3  逐频带事件检测：将连续有源帧拼成事件段
%
%   EventListRaw = detect_band_events_m3(Y_dBm_all, Y_lin_all, Config)
%
%   算法：
%     1. 对每个频带 b、每帧 t，计算总能量 E_b(t) = sum_m Y_lin(m,b,t)
%     2. 计算运行基线 B_b(t)（滑动中值或固定基线）
%     3. 当 E_b(t) - B_b(t) > eta_b 时标记为活跃
%     4. 将相邻活跃帧合并成事件区间 [t_start, t_end]
%     5. 过滤太短的事件（< min_duration）
%
%   输入：
%       Y_dBm_all - (M x B x T) 全时段 RSS 观测 (dBm)
%       Y_lin_all - (M x B x T) 全时段 RSS 观测 (线性 mW)
%       Config    - 配置，需含 Config.m3.detect
%
%   输出：
%       EventListRaw(e) - 原始事件列表，每个元素含：
%         .event_id, .band_id, .time_range, .t_start, .t_end, .duration
%         .obs_segment_dBm (M x L), .obs_segment_lin (M x L)

    [M, B, T] = size(Y_dBm_all);

    % --- 填充默认参数 ---
    det = fill_detect_defaults(Config, B);

    event_id_counter = 0;
    EventListRaw = [];

    for b = 1:B
        % 取该频带全时段线性功率 (M x T)
        Y_lin_b = squeeze(Y_lin_all(:, b, :));  % M x T

        % 1. 计算每帧总能量
        E_b = sum(Y_lin_b, 1);  % 1 x T

        % 2. 计算基线
        baseline = compute_baseline(E_b, det.baseline_mode, det.baseline_window);

        % 3. 阈值检测
        eta = det.eta(b);
        is_active = (E_b - baseline) > eta;

        % 4. 合并连续活跃帧为事件
        segments = merge_active_frames(is_active);

        % 5. 过滤短事件并构建输出
        for k = 1:size(segments, 1)
            ts = segments(k, 1);
            te = segments(k, 2);
            dur = te - ts + 1;

            if dur < det.min_duration
                continue;
            end

            event_id_counter = event_id_counter + 1;
            ev = struct();
            ev.event_id = event_id_counter;
            ev.band_id  = b;
            ev.t_start  = ts;
            ev.t_end    = te;
            ev.time_range    = [ts, te];
            ev.duration      = dur;
            ev.obs_segment_dBm = squeeze(Y_dBm_all(:, b, ts:te));  % M x L
            ev.obs_segment_lin = squeeze(Y_lin_all(:, b, ts:te));   % M x L

            % 确保 M x L 维度（单帧时 squeeze 可能变列向量）
            if dur == 1
                ev.obs_segment_dBm = ev.obs_segment_dBm(:);  % M x 1
                ev.obs_segment_lin = ev.obs_segment_lin(:);
            end

            if isempty(EventListRaw)
                EventListRaw = ev;
            else
                EventListRaw(end+1) = ev; %#ok<AGROW>
            end
        end
    end

    fprintf('[M3] 事件检测完成: %d 个原始事件 (B=%d 频带, T=%d 帧)\n', ...
        numel(EventListRaw), B, T);
end


%% ==================== 局部函数 ====================

function det = fill_detect_defaults(Config, B)
% FILL_DETECT_DEFAULTS  填充 M3 检测参数默认值

    det = struct();

    if isfield(Config, 'm3') && isfield(Config.m3, 'detect')
        d = Config.m3.detect;
    else
        d = struct();
    end

    % 基线模式：'running_median' 或 'fixed'
    if isfield(d, 'baseline_mode')
        det.baseline_mode = d.baseline_mode;
    else
        det.baseline_mode = 'running_median';
    end

    % 基线滑窗长度
    if isfield(d, 'baseline_window')
        det.baseline_window = d.baseline_window;
    else
        det.baseline_window = 20;
    end

    % 每频带检测阈值（线性功率增量）
    if isfield(d, 'eta')
        det.eta = d.eta;
    else
        % 默认阈值：基于典型噪声底估计
        % 噪声功率 ≈ 10^((-174+10*log10(BW))/10)，阈值取噪声 3~5 倍
        det.eta = ones(1, B) * 1e-15;  % 保守默认，可由 Config 覆盖
    end
    if isscalar(det.eta)
        det.eta = repmat(det.eta, 1, B);
    end

    % 最短事件持续帧数
    if isfield(d, 'min_duration')
        det.min_duration = d.min_duration;
    else
        det.min_duration = 1;
    end
end


function baseline = compute_baseline(E_b, mode, window)
% COMPUTE_BASELINE  计算能量基线
%   E_b: 1 x T
%   输出 baseline: 1 x T

    T = numel(E_b);
    baseline = zeros(1, T);

    switch mode
        case 'running_median'
            for t = 1:T
                t_from = max(1, t - window);
                t_to   = t - 1;
                if t_to < t_from
                    baseline(t) = 0;
                else
                    baseline(t) = median(E_b(t_from:t_to));
                end
            end

        case 'fixed'
            % 使用全局中值作为基线
            baseline(:) = median(E_b);

        otherwise
            error('M3:detect', '未知基线模式: %s', mode);
    end
end


function segments = merge_active_frames(is_active)
% MERGE_ACTIVE_FRAMES  将连续 true 帧合并为 [start, end] 区间
%   is_active: 1 x T logical
%   segments: K x 2 矩阵，每行 [t_start, t_end]

    segments = [];
    T = numel(is_active);
    in_event = false;
    t_start = 0;

    for t = 1:T
        if is_active(t) && ~in_event
            in_event = true;
            t_start = t;
        elseif ~is_active(t) && in_event
            in_event = false;
            segments(end+1, :) = [t_start, t-1]; %#ok<AGROW>
        end
    end

    % 处理尾部事件
    if in_event
        segments(end+1, :) = [t_start, T]; %#ok<AGROW>
    end
end
