function GroupList = extract_group_features_m3( ...
    GroupList, EventList, Y_dBm_all, Y_lin_all, SpatialFP, Config)
% EXTRACT_GROUP_FEATURES_M3  group 级特征提取（功率/稳定性/覆盖）
%
%   关键输出：
%     power_excess_mean_dB_per_band
%     power_excess_top1_dB_per_band
%     n_valid_ap_per_band
%     power_zscore_per_band

    %#ok<INUSD>
    if isempty(GroupList)
        return;
    end

    cfg = fill_group_feature_defaults(Config, size(Y_dBm_all, 2));
    B = max(cfg.num_bands, max([EventList.band_id]));
    GroupList = ensure_group_fields_extract(GroupList, B);

    [mu_band, sigma_band] = estimate_band_power_stats(EventList, B, cfg);

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        active_bands = grp.active_bands(:).';

        grp.active_band_count = numel(active_bands);
        grp.band_coverage_ratio = grp.active_band_count / max(B, 1);
        grp.duration_frames = grp.end_frame - grp.start_frame + 1;

        grp.power_excess_mean_dB_per_band = nan(1, B);
        grp.power_excess_top1_dB_per_band = nan(1, B);
        grp.power_stability_per_band = nan(1, B);
        grp.n_valid_ap_per_band = zeros(1, B);
        grp.power_zscore_per_band = nan(1, B);

        band_best_grid_idx = nan(1, B);
        band_best_pos_xy = nan(B, 2);
        noise_floor_per_band = nan(1, B);

        for b = active_bands
            eidx = grp.member_event_idx_per_band(b);
            if eidx <= 0 || eidx > numel(EventList)
                continue;
            end
            ev = EventList(eidx);

            p_excess = get_field_default(ev, 'power_excess_mean_dB', NaN);
            grp.power_excess_mean_dB_per_band(b) = p_excess;
            grp.power_excess_top1_dB_per_band(b) = get_field_default(ev, 'power_excess_top1_dB', NaN);
            grp.power_stability_per_band(b) = get_field_default(ev, 'power_stability_est', NaN);
            grp.n_valid_ap_per_band(b) = get_field_default(ev, 'n_valid_ap', 0);
            noise_floor_per_band(b) = get_field_default(ev, 'noise_floor_dBm', NaN);

            grp.power_zscore_per_band(b) = (p_excess - mu_band(b)) / (sigma_band(b) + cfg.zscore_eps);

            [~, s_obs] = aggregate_event_fingerprint_m4(ev);
            F_shape = SpatialFP.band(b).F_shape_l1;
            dvec = sqrt(sum((s_obs - F_shape).^2, 1));
            [~, g_best] = min(dvec);
            band_best_grid_idx(b) = g_best;
            band_best_pos_xy(b, :) = SpatialFP.grid_xy(g_best, :);
        end

        grp.band_power_mu = mu_band;
        grp.band_power_sigma = sigma_band;
        grp.noise_floor_dBm_per_band = noise_floor_per_band;

        z_active = grp.power_zscore_per_band(active_bands);
        grp.strongest_band_count = sum(z_active >= cfg.strong_band_z);

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
        grp.band_best_pos_xy = band_best_pos_xy;

        GroupList(g) = grp;
    end
end


%% ==================== 局部函数 ====================

function cfg = fill_group_feature_defaults(Config, B)
    cfg.num_bands = B;
    cfg.zscore_eps = 1e-6;
    cfg.min_sigma = 1.0;
    cfg.strong_band_z = 0.8;

    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        cfg.num_bands = Config.m0.num_bands;
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'ordinary') && ...
       isfield(Config.m3.ordinary, 'min_detect_z')
        cfg.strong_band_z = Config.m3.ordinary.min_detect_z;
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


function GroupList = ensure_group_fields_extract(GroupList, B)
    defs = struct();
    defs.band_coverage_ratio = 0;
    defs.power_excess_mean_dB_per_band = nan(1, B);
    defs.power_excess_top1_dB_per_band = nan(1, B);
    defs.power_stability_per_band = nan(1, B);
    defs.n_valid_ap_per_band = zeros(1, B);
    defs.power_zscore_per_band = nan(1, B);
    defs.band_power_mu = zeros(1, B);
    defs.band_power_sigma = ones(1, B);
    defs.noise_floor_dBm_per_band = nan(1, B);
    defs.strongest_band_count = 0;
    defs.multiband_position_spread = NaN;
    defs.multiband_position_consistency = NaN;
    defs.band_best_grid_idx = nan(1, B);
    defs.band_best_pos_xy = nan(B, 2);

    fn = fieldnames(defs);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(GroupList, f)
            [GroupList.(f)] = deal(defs.(f));
        end
    end
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
