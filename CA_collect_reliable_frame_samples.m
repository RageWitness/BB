function CAFrameSamples = CA_collect_reliable_frame_samples( ...
    SourceContext, Y_lin_all, Y_dBm_all, Config)
% CA_COLLECT_RELIABLE_FRAME_SAMPLES  按帧收集 CA 可靠样本
%
%   对每个 (t, b) 检查 label_map 是否为 calibration source，
%   逐帧提取 centered 观测，计算段内/AP/时间权重。

    Config = CA_fill_defaults(Config);
    LC = source_label_constants();
    [M_ap, B, T] = size(Y_lin_all);
    tau_rec = Config.ca.weight.tau_rec;

    label_map = SourceContext.per_band_label_map;   % T x B
    uid_map   = SourceContext.per_band_uid_map;     % T x B cell

    % --- 第 1 遍：收集所有候选帧样本 ---
    CAFrameSamples = [];
    sid = 0;

    for t = 1:T
        for b = 1:B
            lbl = label_map(t, b);
            if lbl ~= LC.PERSISTENT_CAL && lbl ~= LC.BROADBAND_CAL
                continue;
            end

            sid = sid + 1;
            s = struct();
            s.sample_id  = sid;
            s.frame_id   = t;
            s.band_id    = b;
            s.label      = lbl;

            s.source_uid = uid_map{t, b};
            s.location_prior = SourceContext.per_band_location_prior_map{t, b};
            s.power_prior    = SourceContext.per_band_power_prior_map{t, b};
            s.metadata       = SourceContext.per_band_metadata_map{t, b};

            s.valid = true;
            s.reject_reason = '';

            % --- 位置 ---
            if ~isfield(s.location_prior, 'type') || ...
               ~strcmp(s.location_prior.type, 'exact')
                s.valid = false;
                s.reject_reason = 'location_not_exact';
            end

            if s.valid && (isempty(s.location_prior.value))
                s.valid = false;
                s.reject_reason = 'location_value_empty';
            end

            if s.valid
                s.position_xy = s.location_prior.value(:)';
            else
                s.position_xy = [];
            end

            % --- 功率 ---
            if s.valid
                pt = s.power_prior.type;
                if strcmp(pt, 'exact')
                    % ok
                elseif strcmp(pt, 'exact_by_band')
                    if b > numel(s.power_prior.value)
                        s.valid = false;
                        s.reject_reason = 'power_not_exact_for_band';
                    end
                else
                    s.valid = false;
                    s.reject_reason = 'power_not_exact';
                end
            end

            % --- 观测 ---
            y_lin = Y_lin_all(:, b, t);
            if numel(y_lin) ~= M_ap
                s.valid = false;
                s.reject_reason = 'obs_dimension_invalid';
            end

            s.obs_centered_dBm = [];
            s.obs_mean_dBm     = [];
            s.valid_ap_mask    = [];
            s.n_valid_aps      = 0;

            if s.valid
                [z_c, y_dBm, ~, vam, cinfo] = ...
                    CA_compute_centered_observation(y_lin(:), Config);
                s.obs_centered_dBm = z_c;
                s.obs_mean_dBm     = y_dBm;
                s.valid_ap_mask    = vam;
                s.n_valid_aps      = cinfo.n_valid_aps;

                if strcmp(cinfo.status, 'too_few_valid_aps')
                    s.valid = false;
                    s.reject_reason = 'too_few_valid_aps';
                end
            end

            % 权重占位
            s.w_seg = 1;
            s.w_ap  = 0;
            s.w_rec = 0;
            s.sample_weight = 0;

            if s.valid
                s.w_ap  = s.n_valid_aps / M_ap;
                s.w_rec = exp(-(T - t) / tau_rec);
            end

            if isempty(CAFrameSamples)
                CAFrameSamples = s;
            else
                CAFrameSamples(end+1) = s; %#ok<AGROW>
            end
        end
    end

    if isempty(CAFrameSamples)
        CAFrameSamples = struct([]);
        fprintf('[CA-Frame] 无候选帧样本\n');
        return;
    end

    % --- 第 2 遍：计算段内均衡权重 w_seg ---
    valid_mask = [CAFrameSamples.valid];
    valid_idx  = find(valid_mask);

    if isempty(valid_idx)
        fprintf('[CA-Frame] 样本总数=%d  有效=0\n', numel(CAFrameSamples));
        return;
    end

    % group by (source_uid, band_id) → 连续段 → L_seg
    uids   = {CAFrameSamples(valid_idx).source_uid};
    bids   = [CAFrameSamples(valid_idx).band_id];
    fids   = [CAFrameSamples(valid_idx).frame_id];

    seg_keys = cell(1, numel(valid_idx));
    for i = 1:numel(valid_idx)
        seg_keys{i} = sprintf('%s__b%d', uids{i}, bids(i));
    end

    [unique_segs, ~, ic] = unique(seg_keys, 'stable');

    for g = 1:numel(unique_segs)
        g_mask = (ic == g);
        g_frames = sort(fids(g_mask));
        % split on gaps > 1
        breaks = [0, find(diff(g_frames) > 1), numel(g_frames)];
        g_local_idx = find(g_mask);

        for si = 1:numel(breaks)-1
            seg_range = (breaks(si)+1):breaks(si+1);
            L_seg = numel(seg_range);
            for ii = seg_range
                real_idx = valid_idx(g_local_idx(ii));
                CAFrameSamples(real_idx).w_seg = 1 / L_seg;
            end
        end
    end

    % --- 第 3 遍：总权重 ---
    for i = valid_idx
        s = CAFrameSamples(i);
        CAFrameSamples(i).sample_weight = s.w_seg * s.w_ap * s.w_rec;
    end

    n_total = numel(CAFrameSamples);
    n_valid = sum(valid_mask);
    fprintf('[CA-Frame] 帧样本: total=%d  valid=%d\n', n_total, n_valid);
    for b = 1:B
        bm = [CAFrameSamples.band_id] == b & valid_mask;
        if any(bm)
            ws = [CAFrameSamples(bm).sample_weight];
            n_eff_b = sum(ws)^2 / sum(ws.^2);
            fprintf('  band %d: n_raw=%d  n_eff=%.1f  sum_w=%.2f\n', ...
                b, sum(bm), n_eff_b, sum(ws));
        end
    end
end
