function GroupList = extract_group_features_m3( ...
    GroupList, EventList, Y_dBm_all, Y_lin_all, SpatialFP, Config)
% EXTRACT_GROUP_FEATURES_M3  group-level feature extraction (linear-domain primary)

    %#ok<INUSD>
    if isempty(GroupList)
        return;
    end

    [M_data, B_data, T_data] = size(Y_dBm_all);
    cfg = fill_group_feature_defaults(Config, B_data);

    max_event_band = 0;
    if ~isempty(EventList)
        max_event_band = max([EventList.band_id]); %#ok<NBRAK2>
    end
    B = max([cfg.num_bands, B_data, max_event_band]);

    GroupList = ensure_group_fields_extract(GroupList, B, M_data);
    [mu_band, sigma_band] = estimate_band_power_stats(EventList, B, cfg);

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        active_bands = grp.active_bands(:).';

        grp.active_band_count = numel(active_bands);
        grp.band_coverage_ratio = grp.active_band_count / max(B, 1);
        grp.duration_frames = grp.end_frame - grp.start_frame + 1;

        grp.power_excess_mean_dB_per_band = nan(1, B);
        grp.power_excess_top1_dB_per_band = nan(1, B);
        grp.power_stability_per_band      = nan(1, B);
        grp.power_stability_lin_per_band  = nan(1, B);
        grp.n_valid_ap_per_band           = zeros(1, B);
        grp.power_zscore_per_band         = nan(1, B);

        grp.energy_sum_lin_per_band   = zeros(1, B);
        grp.energy_top1_lin_per_band  = zeros(1, B);
        grp.energy_topK_lin_per_band  = zeros(1, B);
        grp.n_valid_ap_lin_per_band   = zeros(1, B);
        grp.noise_floor_lin_per_band  = nan(1, B);
        grp.ratio_sum_lin_per_band    = zeros(1, B);
        grp.ratio_top1_lin_per_band   = zeros(1, B);
        grp.ratio_topK_lin_per_band   = zeros(1, B);

        band_best_grid_idx = nan(1, B);
        band_best_pos_xy   = nan(B, 2);

        for b = active_bands
            eidx = grp.member_event_idx_per_band(b);
            if eidx <= 0 || eidx > numel(EventList)
                continue;
            end
            ev = EventList(eidx);

            p_excess = get_field_default(ev, 'power_excess_mean_dB', NaN);
            grp.power_excess_mean_dB_per_band(b) = p_excess;
            grp.power_excess_top1_dB_per_band(b) = get_field_default(ev, 'power_excess_top1_dB', NaN);
            grp.power_stability_per_band(b)      = get_field_default(ev, 'power_stability_est', NaN);
            grp.n_valid_ap_per_band(b)           = get_field_default(ev, 'n_valid_ap', 0);
            grp.power_zscore_per_band(b)         = (p_excess - mu_band(b)) / (sigma_band(b) + cfg.zscore_eps);

            grp.energy_sum_lin_per_band(b)  = get_field_default(ev, 'energy_sum_lin', 0);
            grp.energy_top1_lin_per_band(b) = get_field_default(ev, 'energy_top1_lin', 0);
            grp.energy_topK_lin_per_band(b) = get_field_default(ev, 'energy_topK_lin', 0);
            grp.n_valid_ap_lin_per_band(b)  = get_field_default(ev, 'n_valid_ap_lin', 0);
            grp.noise_floor_lin_per_band(b) = get_field_default(ev, 'noise_floor_lin', NaN);
            grp.ratio_sum_lin_per_band(b)   = get_field_default(ev, 'ratio_sum_lin', 0);
            grp.ratio_top1_lin_per_band(b)  = get_field_default(ev, 'ratio_top1_lin', 0);
            grp.ratio_topK_lin_per_band(b)  = get_field_default(ev, 'ratio_topK_lin', 0);
            grp.power_stability_lin_per_band(b) = get_field_default(ev, 'power_stability_lin', 0);

            [~, s_obs] = aggregate_event_fingerprint_m4(ev);
            F_shape = SpatialFP.band(b).F_shape_l1;
            dvec = sqrt(sum((s_obs - F_shape).^2, 1));
            [~, g_best] = min(dvec);
            band_best_grid_idx(b) = g_best;
            band_best_pos_xy(b, :) = SpatialFP.grid_xy(g_best, :);
        end

        grp.band_power_mu = mu_band;
        grp.band_power_sigma = sigma_band;

        valid_pos = ~isnan(band_best_grid_idx(active_bands));
        if sum(valid_pos) >= 2
            pos_use = band_best_pos_xy(active_bands(valid_pos), :);
            c = mean(pos_use, 1);
            d = sqrt(sum((pos_use - c).^2, 2));
            spread = max(d);
        else
            spread = NaN;
        end
        grp.multiband_position_spread = spread;
        grp.multiband_position_consistency = spread;
        grp.band_best_grid_idx = band_best_grid_idx;
        grp.band_best_pos_xy   = band_best_pos_xy;

        z_active = grp.power_zscore_per_band(active_bands);
        grp.strongest_band_count = sum(z_active(~isnan(z_active)) >= cfg.strong_band_z);

        [pwr_mat, min_pwr] = build_group_power_matrix_lin( ...
            Y_lin_all, grp.start_frame, grp.end_frame, B, B_data, T_data, cfg.trusted_hard.aggregation);
        grp.trusted_hard_power_matrix_lin = pwr_mat;
        grp.trusted_hard_min_power_lin = min_pwr;
        grp.trusted_hard_active_band_count = grp.active_band_count;
        grp.trusted_hard_duration_frames = grp.duration_frames;
        grp.trusted_hard_pass = false;
        grp.trusted_hard_removed_from_scoring = false;
        grp.trusted_hard_aggregation = cfg.trusted_hard.aggregation;

        GroupList(g) = grp;
    end
end


%% ==================== local helpers ====================

function cfg = fill_group_feature_defaults(Config, B)
    cfg.num_bands = B;
    cfg.zscore_eps = 1e-6;
    cfg.min_sigma = 1.0;
    cfg.strong_band_z = 0.8;
    cfg.trusted_hard.aggregation = 'mean';

    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        cfg.num_bands = Config.m0.num_bands;
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'ordinary') && ...
       isfield(Config.m3.ordinary, 'min_detect_z')
        cfg.strong_band_z = Config.m3.ordinary.min_detect_z;
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'trusted_hard')
        th = Config.m3.trusted_hard;
        if isfield(th, 'aggregation') && ischar(th.aggregation)
            cfg.trusted_hard.aggregation = lower(strtrim(th.aggregation));
        end
    end
end


function [mu_band, sigma_band] = estimate_band_power_stats(EventList, B, cfg)
    mu_band = zeros(1, B);
    sigma_band = ones(1, B) * cfg.min_sigma;

    for b = 1:B
        idx = find([EventList.band_id] == b); %#ok<NBRAK2>
        if isempty(idx)
            continue;
        end
        v = arrayfun(@(e)get_field_default(e, 'power_excess_mean_dB', NaN), EventList(idx));
        v = v(~isnan(v) & isfinite(v));
        if isempty(v)
            continue;
        end
        mu_band(b) = mean(v);
        s = std(v);
        if ~isfinite(s) || s < cfg.min_sigma
            s = cfg.min_sigma;
        end
        sigma_band(b) = s;
    end
end


function [pwr_mat, min_pwr] = build_group_power_matrix_lin( ...
    Y_lin_all, t_start, t_end, B, B_data, T_data, aggregation)

    M = size(Y_lin_all, 1);
    pwr_mat = nan(M, B);
    min_pwr = NaN;

    if isempty(Y_lin_all) || T_data <= 0
        return;
    end

    ts = max(1, round(t_start));
    te = min(T_data, round(t_end));
    if te < ts
        return;
    end

    b_use = min(B, B_data);
    for b = 1:b_use
        seg = squeeze(Y_lin_all(:, b, ts:te));  % M x L or M x 1
        if isempty(seg)
            continue;
        end
        if isvector(seg)
            seg = seg(:);
        end
        switch aggregation
            case 'median'
                pwr_mat(:, b) = median(seg, 2);
            case 'min'
                pwr_mat(:, b) = min(seg, [], 2);
            otherwise
                pwr_mat(:, b) = mean(seg, 2); % default: mean
        end
    end

    min_pwr = nanmin_local(pwr_mat(:));
end


function GroupList = ensure_group_fields_extract(GroupList, B, M)
    defs = struct();
    defs.band_coverage_ratio = 0;
    defs.power_excess_mean_dB_per_band = nan(1, B);
    defs.power_excess_top1_dB_per_band = nan(1, B);
    defs.power_stability_per_band      = nan(1, B);
    defs.power_stability_lin_per_band  = nan(1, B);
    defs.n_valid_ap_per_band           = zeros(1, B);
    defs.power_zscore_per_band         = nan(1, B);
    defs.band_power_mu    = zeros(1, B);
    defs.band_power_sigma = ones(1, B);
    defs.strongest_band_count = 0;
    defs.multiband_position_spread = NaN;
    defs.multiband_position_consistency = NaN;
    defs.band_best_grid_idx = nan(1, B);
    defs.band_best_pos_xy = nan(B, 2);
    defs.energy_sum_lin_per_band = zeros(1, B);
    defs.energy_top1_lin_per_band = zeros(1, B);
    defs.energy_topK_lin_per_band = zeros(1, B);
    defs.n_valid_ap_lin_per_band = zeros(1, B);
    defs.noise_floor_lin_per_band = nan(1, B);
    defs.ratio_sum_lin_per_band = zeros(1, B);
    defs.ratio_top1_lin_per_band = zeros(1, B);
    defs.ratio_topK_lin_per_band = zeros(1, B);
    defs.trusted_hard_power_matrix_lin = nan(M, B);
    defs.trusted_hard_min_power_lin = NaN;
    defs.trusted_hard_active_band_count = 0;
    defs.trusted_hard_duration_frames = 0;
    defs.trusted_hard_pass = false;
    defs.trusted_hard_removed_from_scoring = false;
    defs.trusted_hard_aggregation = 'mean';

    fn = fieldnames(defs);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(GroupList, f)
            [GroupList.(f)] = deal(defs.(f));
        end
    end
end


function v = nanmin_local(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = min(x);
    end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
