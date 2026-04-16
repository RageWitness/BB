function EventList = extract_event_features_m3(EventListRaw, Y_dBm_all, Y_lin_all, ...
    SpatialFP, SignatureLib, SourceTemplates, Config)
% EXTRACT_EVENT_FEATURES_M3  对每个原始事件提取分类特征
%
%   保留原有特征字段，同时新增相对噪声底功率特征：
%     power_excess_mean_dB
%     power_excess_top1_dB
%     power_excess_topk_dB
%     n_valid_ap

    %#ok<INUSD>
    [~, B, ~] = size(Y_dBm_all);
    n_events = numel(EventListRaw);

    % 填充默认参数
    feat_cfg = fill_feat_defaults(Config);

    % 线性域配置
    lp_cfg = fill_linear_power_defaults(Config);

    % 先为所有事件预分配字段（MATLAB struct 数组要求字段一致）
    for e = 1:n_events
        EventListRaw(e).power_level_est            = 0;
        EventListRaw(e).power_stability_est        = 0;
        EventListRaw(e).power_stability_lin        = 0;
        EventListRaw(e).band_coverage_vec          = zeros(1, B);
        EventListRaw(e).schedule_match_score       = 0;
        EventListRaw(e).position_prior_match_score = 0;
        EventListRaw(e).power_excess_mean_dB       = 0;
        EventListRaw(e).power_excess_top1_dB       = 0;
        EventListRaw(e).power_excess_topk_dB       = 0;
        EventListRaw(e).n_valid_ap                 = 0;
        EventListRaw(e).noise_floor_dBm            = 0;
        EventListRaw(e).energy_sum_lin             = 0;
        EventListRaw(e).energy_top1_lin            = 0;
        EventListRaw(e).energy_topK_lin            = 0;
        EventListRaw(e).n_valid_ap_lin             = 0;
        EventListRaw(e).noise_floor_lin            = 0;
        EventListRaw(e).ratio_sum_lin              = 0;
        EventListRaw(e).ratio_top1_lin             = 0;
        EventListRaw(e).ratio_topK_lin             = 0;
        EventListRaw(e).feat_vec                   = [];
    end
    EventList = EventListRaw;

    for e = 1:n_events
        ev = EventList(e);
        b  = ev.band_id;
        L  = ev.duration;

        %% 1) 功率等级估计（旧字段保留）
        % 每帧 AP 平均 RSS (dBm)
        P_bar_per_frame = mean(ev.obs_segment_dBm, 1);  % 1 x L
        ev.power_level_est = mean(P_bar_per_frame);

        %% 2) 功率稳定性估计（旧字段保留）
        if L > 1
            ev.power_stability_est = var(P_bar_per_frame);
        else
            ev.power_stability_est = 0;
        end

        %% 3) 频带覆盖向量
        ev.band_coverage_vec = compute_band_coverage(ev, EventListRaw, B);

        %% 4) 时间先验匹配分数
        ev.schedule_match_score = compute_schedule_match( ...
            ev, SignatureLib, SourceTemplates, Config);

        %% 5) 位置先验快速匹配分数
        ev.position_prior_match_score = compute_position_match( ...
            ev, SpatialFP, SignatureLib, feat_cfg);

        %% 6) 相对噪声底功率特征（新增）
        rss_mean_per_ap_dBm = mean(ev.obs_segment_dBm, 2);  % M x 1
        noise_floor_dBm = get_noise_floor_dBm_for_band_m3( ...
            Config, b, numel(rss_mean_per_ap_dBm));
        excess_dB = rss_mean_per_ap_dBm - noise_floor_dBm;

        valid_mask = excess_dB >= feat_cfg.tau_valid_ap_dB;
        ev.n_valid_ap = sum(valid_mask);
        ev.noise_floor_dBm = noise_floor_dBm(1);

        if any(valid_mask)
            ev.power_excess_mean_dB = mean(excess_dB(valid_mask));
        else
            ev.power_excess_mean_dB = mean(excess_dB);
        end

        ev.power_excess_top1_dB = max(excess_dB);
        k_use = min(feat_cfg.topk, numel(excess_dB));
        ex_sorted = sort(excess_dB, 'descend');
        ev.power_excess_topk_dB = mean(ex_sorted(1:k_use));

        %% 7) 线性域功率特征（新增）
        if ev.duration > 1
            y_lin_mean = mean(ev.obs_segment_lin, 2);  % M x 1
        else
            y_lin_mean = ev.obs_segment_lin(:);
        end
        nf_lin = 10^(noise_floor_dBm(1) / 10);  % 标量噪声底（线性 mW）
        nf_lin = max(nf_lin, 1e-30);
        ev.noise_floor_lin = nf_lin;

        % 有效 AP：线性功率 > alpha_noise * noise_floor
        valid_lin_mask = y_lin_mean > (lp_cfg.alpha_noise * nf_lin);
        ev.n_valid_ap_lin = sum(valid_lin_mask);

        y_valid = y_lin_mean(valid_lin_mask);
        ev.energy_sum_lin = sum(y_valid);
        ev.energy_top1_lin = max(y_lin_mean);

        K_use = min(lp_cfg.topk, numel(y_lin_mean));
        y_sorted = sort(y_lin_mean, 'descend');
        ev.energy_topK_lin = sum(y_sorted(1:K_use));

        N_eff = max(ev.n_valid_ap_lin, 1);
        ev.ratio_sum_lin  = ev.energy_sum_lin / (nf_lin * N_eff);
        ev.ratio_top1_lin = ev.energy_top1_lin / nf_lin;
        ev.ratio_topK_lin = ev.energy_topK_lin / (nf_lin * K_use);

        % 线性域功率稳定性（帧间方差 / 帧间均值）
        if ev.duration > 1
            frame_energy = sum(ev.obs_segment_lin, 1);  % 1 x L
            ev.power_stability_lin = var(frame_energy) / (mean(frame_energy)^2 + 1e-30);
        else
            ev.power_stability_lin = 0;
        end

        %% 8) 组合特征向量（旧 + 新）
        ev.feat_vec = [ev.power_level_est, ev.power_stability_est, L, ...
                       ev.band_coverage_vec, ...
                       ev.schedule_match_score, ev.position_prior_match_score, ...
                       ev.power_excess_mean_dB, ev.power_excess_top1_dB, ...
                       ev.power_excess_topk_dB, ev.n_valid_ap, ...
                       ev.ratio_sum_lin, ev.ratio_top1_lin, ev.ratio_topK_lin, ...
                       ev.n_valid_ap_lin];

        EventList(e) = ev;
    end

    fprintf('[M3] 事件特征提取完成: %d 个事件\n', n_events);
end


%% ==================== 局部函数 ====================

function cfg = fill_feat_defaults(Config)
% FILL_FEAT_DEFAULTS  填充特征提取参数默认值
    cfg = struct();
    if isfield(Config, 'm3') && isfield(Config.m3, 'feat')
        f = Config.m3.feat;
    else
        f = struct();
    end

    % 位置匹配高斯核带宽
    if isfield(f, 'sigma_pos')
        cfg.sigma_pos = f.sigma_pos;
    else
        cfg.sigma_pos = 10;  % dB
    end

    if isfield(Config, 'm3') && isfield(Config.m3, 'relative_power')
        rp = Config.m3.relative_power;
    else
        rp = struct();
    end

    if isfield(rp, 'topk')
        cfg.topk = rp.topk;
    else
        cfg.topk = 3;
    end
    if isfield(rp, 'tau_valid_ap_dB')
        cfg.tau_valid_ap_dB = rp.tau_valid_ap_dB;
    else
        cfg.tau_valid_ap_dB = 3;
    end
end


function coverage = compute_band_coverage(ev, EventListRaw, B)
% COMPUTE_BAND_COVERAGE  计算跨频带同步覆盖向量
%   如果其他频带上有与当前事件时间重叠的事件，标记为 1
%
%   coverage: 1 x B, coverage(b)=1 表示该频带有同步事件

    coverage = zeros(1, B);
    coverage(ev.band_id) = 1;  % 自身频带

    for k = 1:numel(EventListRaw)
        other = EventListRaw(k);
        if other.band_id == ev.band_id
            continue;
        end
        % 检查时间重叠
        overlap_start = max(ev.t_start, other.t_start);
        overlap_end   = min(ev.t_end,   other.t_end);
        if overlap_end >= overlap_start
            coverage(other.band_id) = 1;
        end
    end
end


function score = compute_schedule_match(ev, SignatureLib, SourceTemplates, Config)
% COMPUTE_SCHEDULE_MATCH  计算事件与先验模板的最佳时间匹配分数
%   score = max_j S^time_j(e) = |T_e ∩ T_j^prior| / (|T_e| + eps)

    score = 0;
    n_tpl = SignatureLib.n_templates;
    event_frames = ev.t_start : ev.t_end;
    L = ev.duration;

    for j = 1:n_tpl
        tpl = SignatureLib.templates(j);

        % 只匹配有时间先验的模板（scheduled 类型）
        if ~strcmp(tpl.time_pattern_type, 'scheduled')
            continue;
        end

        % 获取完整时间先验信息
        time_prior = get_time_prior_from_templates(tpl, SourceTemplates, Config);
        if isempty(time_prior)
            continue;
        end

        % 计算预期活跃帧集合
        expected_frames = get_expected_active_frames(time_prior, Config.m0.T_total);

        % 交集占比
        overlap = numel(intersect(event_frames, expected_frames));
        s_j = overlap / (L + 1e-10);
        score = max(score, s_j);
    end
end


function time_prior = get_time_prior_from_templates(tpl_entry, SourceTemplates, Config)
% GET_TIME_PRIOR_FROM_TEMPLATES  根据 SignatureLib 条目反查原始模板的 time_prior
    time_prior = [];
    key = tpl_entry.template_key;

    % 查 prior_pos
    B = Config.m0.num_bands;
    for b = 1:B
        Sprior_b = Config.m0.prior.Sprior(b);
        for s = 1:Sprior_b
            if strcmp(SourceTemplates.prior_pos(b,s).template_key, key)
                time_prior = SourceTemplates.prior_pos(b,s).time_prior;
                return;
            end
        end
    end

    % 查 prior_time
    for b = 1:B
        Tprior_b = Config.m0.prior.Tprior(b);
        for q = 1:Tprior_b
            if strcmp(SourceTemplates.prior_time(b,q).template_key, key)
                time_prior = SourceTemplates.prior_time(b,q).time_prior;
                return;
            end
        end
    end
end


function frames = get_expected_active_frames(time_prior, T_total)
% GET_EXPECTED_ACTIVE_FRAMES  根据 time_prior 生成预期活跃帧集合
    frames = [];

    switch time_prior.mode
        case 'periodic'
            period = time_prior.period_frames;
            dur    = time_prior.duration_frames;
            phase  = time_prior.phase_frames;
            for t = 1:T_total
                p = mod(t - 1 - phase, period);
                if p < dur
                    frames(end+1) = t; %#ok<AGROW>
                end
            end

        case 'table'
            if isfield(time_prior, 'schedule') && ~isempty(time_prior.schedule)
                frames = time_prior.schedule;
            end

        otherwise
            % unknown 模式，无预期帧
    end
end


function score = compute_position_match(ev, SpatialFP, SignatureLib, feat_cfg)
% COMPUTE_POSITION_MATCH  计算事件与位置已知模板的最佳指纹匹配分数
%   对位置已知模板，用功率修正后的快速距离计算相似度

    score = 0;
    b = ev.band_id;
    n_tpl = SignatureLib.n_templates;

    % 事件平均观测指纹 (M x 1)
    if ev.duration > 1
        F_obs = mean(ev.obs_segment_dBm, 2);  % M x 1
    else
        F_obs = ev.obs_segment_dBm(:);
    end

    % 观测中心化
    mu_obs = mean(F_obs);
    F_obs_centered = F_obs - mu_obs;

    sigma_pos = feat_cfg.sigma_pos;

    for j = 1:n_tpl
        tpl = SignatureLib.templates(j);

        % 只匹配位置已知的模板
        if ~strcmp(tpl.known_position_mode, 'fixed_known')
            continue;
        end
        if isempty(tpl.candidate_positions)
            continue;
        end

        % 候选位置在指纹库中的最近网格点
        cand_pos = tpl.candidate_positions;  % K x 2
        n_cand = size(cand_pos, 1);

        d_min = inf;
        for c = 1:n_cand
            % 找最近网格点
            dists_to_grid = sqrt(sum((SpatialFP.grid_xy - cand_pos(c,:)).^2, 2));
            [~, g_idx] = min(dists_to_grid);

            % 参考指纹中心化向量
            F_db_centered = SpatialFP.band(b).centered_dBm(:, g_idx);

            % 欧氏距离
            d = norm(F_obs_centered - F_db_centered, 2);
            d_min = min(d_min, d);
        end

        % 高斯核相似度
        s_j = exp(-(d_min^2) / (sigma_pos^2));
        score = max(score, s_j);
    end
end


function noise_floor_dBm = get_noise_floor_dBm_for_band_m3(Config, band_id, M)
% GET_NOISE_FLOOR_DBM_FOR_BAND_M3  M3 特征提取用噪声底
    if nargin < 3 || isempty(M)
        M = 1;
    end

    % 模式一：直接给每个 band 的噪声功率
    if isfield(Config, 'm1') && isfield(Config.m1, 'noise') && ...
       isfield(Config.m1.noise, 'mode') && strcmp(Config.m1.noise.mode, 'power_dBm')
        if isfield(Config.m1.noise, 'noise_power_dBm') && ...
           numel(Config.m1.noise.noise_power_dBm) >= band_id
            nf = Config.m1.noise.noise_power_dBm(band_id);
        else
            nf = -174;
        end
        noise_floor_dBm = repmat(nf, M, 1);
        return;
    end

    % 模式二：n0_dBmHz + BW
    if isfield(Config, 'm1') && isfield(Config.m1, 'noise') && ...
       isfield(Config.m1.noise, 'n0_dBmHz')
        n0_arr = Config.m1.noise.n0_dBmHz;
        if isscalar(n0_arr)
            n0 = n0_arr;
        elseif numel(n0_arr) >= band_id
            n0 = n0_arr(band_id);
        else
            n0 = -174;
        end
    else
        n0 = -174;
    end

    if isfield(Config, 'm1') && isfield(Config.m1, 'bands') && ...
       isfield(Config.m1.bands, 'bw_Hz') && numel(Config.m1.bands.bw_Hz) >= band_id
        bw_hz = Config.m1.bands.bw_Hz(band_id);
    else
        bw_hz = 1;
    end

    nf = n0 + 10 * log10(bw_hz);
    noise_floor_dBm = repmat(nf, M, 1);
end


function cfg = fill_linear_power_defaults(Config)
% FILL_LINEAR_POWER_DEFAULTS  线性域功率特征参数
    cfg.topk = 3;
    cfg.alpha_noise = 2;

    if isfield(Config, 'm3') && isfield(Config.m3, 'linear_power')
        lp = Config.m3.linear_power;
        if isfield(lp, 'topk'),       cfg.topk = lp.topk; end
        if isfield(lp, 'alpha_noise'), cfg.alpha_noise = lp.alpha_noise; end
    end
end
