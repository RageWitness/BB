function [PendingSpatialFP, CAResult] = CA_calibrate_spatialfp_trusted_sources( ...
    SpatialFP, SourceContext, Y_dBm_all, APs, Config)
% CA_CALIBRATE_SPATIALFP_TRUSTED_SOURCES  Build pending RF_raw map from trusted sources.
%
% The input SpatialFP is not modified. The returned PendingSpatialFP is the
% calibrated candidate map for later localization or manual acceptance.

    if nargin < 5
        Config = struct();
    end

    cfg = CA_trusted_rss_defaults(Config);
    LC = source_label_constants();

    PendingSpatialFP = SpatialFP;
    PendingSpatialFP.is_pending_ca_map = true;
    PendingSpatialFP.ca_method = 'trusted_rss_lognormal_top2';

    CAResult = struct();
    CAResult.status = 'ok';
    CAResult.config = cfg;
    CAResult.n_labels_total = 0;
    CAResult.n_candidates_total = 0;
    CAResult.n_candidates_valid = 0;
    CAResult.band = init_band_result(SpatialFP.B);
    CAResult.candidates = struct([]);

    if ~isfield(SourceContext, 'label_table') || isempty(SourceContext.label_table)
        CAResult.status = 'no_labels';
        return;
    end

    [M_obs, B_obs, T_obs] = size(Y_dBm_all);
    if M_obs ~= APs.num
        error('[CA] Y_dBm_all AP dimension does not match APs.num');
    end

    CAResult.n_labels_total = numel(SourceContext.label_table);

    candidate_maps = cell(1, SpatialFP.B);
    candidate_info = cell(1, SpatialFP.B);

    for idx = 1:numel(SourceContext.label_table)
        lbl = SourceContext.label_table(idx);

        if lbl.label ~= LC.PERSISTENT_CAL && lbl.label ~= LC.BROADBAND_CAL
            continue;
        end

        b = lbl.band_id;
        if b < 1 || b > SpatialFP.B || b > B_obs
            cand = rejected_candidate(lbl, 'band_out_of_range');
            CAResult = append_candidate(CAResult, cand);
            continue;
        end

        [ok_pos, source_xy] = get_exact_position(lbl);
        if ~ok_pos
            cand = rejected_candidate(lbl, 'location_not_exact');
            CAResult = append_candidate(CAResult, cand);
            continue;
        end

        sf = max(1, round(lbl.start_frame));
        ef = min(T_obs, round(lbl.end_frame));
        if ef < sf || (ef - sf + 1) < cfg.min_frames
            cand = rejected_candidate(lbl, 'too_few_frames');
            CAResult = append_candidate(CAResult, cand);
            continue;
        end

        rss = squeeze(Y_dBm_all(:, b, sf:ef));  % M x T before reshape
        if isvector(rss)
            rss = rss(:)';                       % 1 x M for a single frame
        elseif size(rss, 1) == M_obs
            rss = rss';                          % T x M
        end

        F_cur = get_rf_raw_map(SpatialFP.band(b));
        P_off = get_offline_power(SpatialFP, b);
        [power_mode, P_c, p_ok, p_reason] = resolve_power_mode(lbl, b, P_off, cfg);
        if ~p_ok
            cand = rejected_candidate(lbl, p_reason);
            CAResult = append_candidate(CAResult, cand);
            continue;
        end

        cfg_one = cfg;
        cfg_one.P_off = P_off;
        if isfinite(P_c)
            cfg_one.P_c = P_c;
        end

        cal = CA_calibrate_fp_by_trusted_source( ...
            APs.pos_xy, SpatialFP.grid_xy, source_xy, rss, F_cur, power_mode, cfg_one);

        cand = struct();
        cand.valid = cal.valid;
        cand.reject_reason = '';
        cand.source_uid = lbl.source_uid;
        cand.band_id = b;
        cand.label = lbl.label;
        cand.frame_range = [sf, ef];
        cand.source_xy = source_xy;
        cand.power_mode = power_mode;
        cand.P_c = P_c;
        cand.P_off = P_off;
        cand.n_hat = cal.n_hat;
        cand.delta_pow = cal.delta_pow;
        cand.gamma_hat = cal.gamma_hat;
        cand.warning = cal.warning;
        cand.search_result = cal.search_result;

        if cal.valid
            if isempty(candidate_maps{b})
                candidate_maps{b} = {cal.F_new};
                candidate_info{b} = {cand};
            else
                candidate_maps{b}{end+1} = cal.F_new; %#ok<AGROW>
                candidate_info{b}{end+1} = cand; %#ok<AGROW>
            end
        else
            cand.reject_reason = cal.warning;
        end

        CAResult = append_candidate(CAResult, cand);
    end

    for b = 1:SpatialFP.B
        maps_b = candidate_maps{b};
        if isempty(maps_b)
            CAResult.band(b).status = 'no_valid_trusted_source';
            continue;
        end

        F_new = fuse_candidate_maps(maps_b, cfg.fusion_mode);
        PendingSpatialFP.band(b).RF_raw = F_new;
        PendingSpatialFP.band(b).F_dBm = F_new;
        PendingSpatialFP.band(b).F_lin = 10.^(F_new / 10);
        PendingSpatialFP.band(b).ca_trusted_rss_applied = true;
        PendingSpatialFP.band(b).ca_trusted_rss_n_candidates = numel(maps_b);
        PendingSpatialFP.band(b).ca_trusted_rss_info = [candidate_info{b}{:}];

        CAResult.band(b).status = 'ok';
        CAResult.band(b).n_candidates = numel(maps_b);
        CAResult.band(b).mean_delta_dB = mean_abs_finite(F_new - get_rf_raw_map(SpatialFP.band(b)));
        CAResult.band(b).n_hat_values = collect_n_hat(candidate_info{b});
    end

    CAResult.n_candidates_total = numel(CAResult.candidates);
    if CAResult.n_candidates_total > 0
        CAResult.n_candidates_valid = sum([CAResult.candidates.valid]);
    end

    if CAResult.n_candidates_valid == 0
        CAResult.status = 'no_valid_trusted_source';
    end

    if cfg.verbose
        fprintf('[CA] trusted RF_raw calibration: valid %d / %d candidates\n', ...
            CAResult.n_candidates_valid, CAResult.n_candidates_total);
        for b = 1:SpatialFP.B
            fprintf('  band %d: %s, candidates=%d\n', ...
                b, CAResult.band(b).status, CAResult.band(b).n_candidates);
        end
    end
end


function band = init_band_result(B)
    band = struct('status', {}, 'n_candidates', {}, 'mean_delta_dB', {}, 'n_hat_values', {});
    for b = 1:B
        band(b).status = 'not_processed';
        band(b).n_candidates = 0;
        band(b).mean_delta_dB = NaN;
        band(b).n_hat_values = [];
    end
end


function [ok, xy] = get_exact_position(lbl)
    ok = false;
    xy = [];
    if isfield(lbl, 'location_prior') && isstruct(lbl.location_prior) && ...
            isfield(lbl.location_prior, 'type') && strcmp(lbl.location_prior.type, 'exact') && ...
            isfield(lbl.location_prior, 'value') && numel(lbl.location_prior.value) == 2
        xy = lbl.location_prior.value(:)';
        ok = all(isfinite(xy));
        return;
    end
    if isfield(lbl, 'position_hint') && numel(lbl.position_hint) == 2
        xy = lbl.position_hint(:)';
        ok = all(isfinite(xy));
    end
end


function [mode, P_c, ok, reason] = resolve_power_mode(lbl, band_id, P_off, cfg)
    mode = 'unknown';
    P_c = NaN;
    ok = true;
    reason = '';

    if ~isfield(lbl, 'power_prior') || ~isstruct(lbl.power_prior) || ...
            ~isfield(lbl.power_prior, 'type')
        return;
    end

    pp = lbl.power_prior;
    switch pp.type
        case 'none'
            mode = 'unknown';
        case 'exact'
            if ~isfield(pp, 'value') || isempty(pp.value)
                ok = false;
                reason = 'power_exact_value_empty';
                return;
            end
            if isscalar(pp.value)
                P_c = pp.value;
            elseif band_id <= numel(pp.value)
                P_c = pp.value(band_id);
            else
                ok = false;
                reason = 'power_exact_value_missing_band';
                return;
            end
        case 'exact_by_band'
            if ~isfield(pp, 'value') || band_id > numel(pp.value)
                ok = false;
                reason = 'power_not_exact_for_band';
                return;
            end
            P_c = pp.value(band_id);
        otherwise
            ok = false;
            reason = 'unsupported_power_prior';
            return;
    end

    if ~isnan(P_c) && (~isscalar(P_c) || ~isfinite(P_c))
        ok = false;
        reason = 'invalid_power_value';
        return;
    end

    if isfinite(P_c)
        if abs(P_c - P_off) <= cfg.power_equal_tol_dB
            mode = 'known_same';
        else
            mode = 'known_different';
        end
    end
end


function P_off = get_offline_power(SpatialFP, band_id)
    if isfield(SpatialFP.band(band_id), 'reference_power_dBm_used')
        P_off = SpatialFP.band(band_id).reference_power_dBm_used;
    elseif isfield(SpatialFP, 'reference_power_dBm_used')
        P = SpatialFP.reference_power_dBm_used;
        P_off = P(min(band_id, numel(P)));
    elseif isfield(SpatialFP, 'ref_power_dBm')
        P = SpatialFP.ref_power_dBm;
        P_off = P(min(band_id, numel(P)));
    else
        error('[CA] Missing offline reference power for band %d', band_id);
    end

    if ~isscalar(P_off) || ~isfinite(P_off)
        error('[CA] Invalid offline reference power for band %d', band_id);
    end
end


function F = get_rf_raw_map(band_fp)
    if isfield(band_fp, 'RF_raw')
        F = band_fp.RF_raw;
    elseif isfield(band_fp, 'F_dBm')
        F = band_fp.F_dBm;
    else
        error('[CA] SpatialFP.band is missing RF_raw/F_dBm');
    end
end


function F = fuse_candidate_maps(maps, mode)
    M = size(maps{1}, 1);
    G = size(maps{1}, 2);
    K = numel(maps);
    X = nan(M, G, K);
    for k = 1:K
        X(:, :, k) = maps{k};
    end

    switch lower(mode)
        case 'median'
            F = median(X, 3, 'omitnan');
        case 'mean'
            num = sum(replace_nan(X, 0), 3);
            den = sum(isfinite(X), 3);
            F = num ./ max(den, 1);
            F(den == 0) = NaN;
        otherwise
            error('[CA] Unknown fusion_mode=%s', mode);
    end
end


function X = replace_nan(X, val)
    X(~isfinite(X)) = NaN;
    X(isnan(X)) = val;
end


function v = mean_abs_finite(X)
    x = X(:);
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = mean(abs(x));
    end
end


function arr = collect_n_hat(info_cell)
    arr = [];
    for k = 1:numel(info_cell)
        if isfield(info_cell{k}, 'n_hat') && isfinite(info_cell{k}.n_hat)
            arr(end+1) = info_cell{k}.n_hat; %#ok<AGROW>
        end
    end
end


function cand = rejected_candidate(lbl, reason)
    cand = struct();
    cand.valid = false;
    cand.reject_reason = reason;
    cand.source_uid = lbl.source_uid;
    cand.band_id = lbl.band_id;
    cand.label = lbl.label;
    cand.frame_range = [lbl.start_frame, lbl.end_frame];
    cand.source_xy = [];
    cand.power_mode = '';
    cand.P_c = NaN;
    cand.P_off = NaN;
    cand.n_hat = NaN;
    cand.delta_pow = NaN;
    cand.gamma_hat = NaN;
    cand.warning = reason;
    cand.search_result = struct();
end


function CAResult = append_candidate(CAResult, cand)
    if isempty(CAResult.candidates)
        CAResult.candidates = cand;
    else
        CAResult.candidates(end+1) = cand;
    end
end
