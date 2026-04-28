function [PendingSpatialFP, CDResult] = CD_calibrate_spatialfp_opportunistic_uncertain( ...
    SpatialFP, calibration_events, APs, Config)
% CD_CALIBRATE_SPATIALFP_OPPORTUNISTIC_UNCERTAIN  Calibrate RF_raw using uncertain opportunistic sources.
%
% This module supports power-known opportunistic events whose location prior
% is a rectangle region or a Gaussian distribution. For each band it first
% estimates n0 candidates from the nearest AP of each usable event, takes a
% band-level median n0, then reuses that n0 for all AP profiles in that band.

    if nargin < 4
        Config = struct();
    end
    cfg = CD_uncertain_rss_defaults(Config);

    PendingSpatialFP = SpatialFP;
    PendingSpatialFP.is_pending_cd_map = true;
    PendingSpatialFP.cd_method = 'opportunistic_uncertain_ray_lognormal';

    CDResult = init_result(SpatialFP);
    CDResult.config = cfg;
    CDResult.n_events_total = numel(calibration_events);

    if isempty(calibration_events)
        CDResult.status = 'no_events';
        return;
    end

    B = SpatialFP.B;
    M = SpatialFP.M;
    G = SpatialFP.G;
    grid_xy = SpatialFP.grid_xy;

    num = cell(1, B);
    den = cell(1, B);
    F_old_cache = cell(1, B);
    for b = 1:B
        F_old_cache{b} = get_band_dbm(SpatialFP.band(b));
        num{b} = zeros(M, G);
        den{b} = cfg.blend_old_weight * ones(M, G);
    end

    [n0_by_band, n0_values] = estimate_shared_n0_by_band( ...
        calibration_events, SpatialFP, F_old_cache, APs, grid_xy, cfg);
    for b = 1:B
        CDResult.band(b).n0_hat_values = n0_values{b};
        CDResult.band(b).n0_hat_median = n0_by_band(b);
    end

    for e = 1:numel(calibration_events)
        ev = calibration_events(e);
        [ok_ev, reason] = accept_event(ev, SpatialFP, cfg);
        if ~ok_ev
            CDResult = append_candidate(CDResult, rejected_candidate(ev, reason));
            continue;
        end

        b = ev.band_id;
        F_old = F_old_cache{b};
        P0 = get_offline_power(SpatialFP, b);
        P_e = resolve_power(ev, b);
        y_tilde = ev.rss_dBm_by_ap(:)' - (P_e - P0);
        q = CB_make_location_distribution_from_prior(ev.location_prior);

        cand = candidate_base(ev);
        cand.P0_dBm = P0;
        cand.P_e_dBm = P_e;
        cand.location_prior_type = q.type;
        cand.n0_mode = cfg.n0_mode;
        n_band_hat = n0_by_band(b);
        cand.n_band_hat = n_band_hat;
        cand.n0_ref_ap = nearest_ap_for_prior(q, APs);
        ref_ap = cand.n0_ref_ap;
        if ~isfinite(n_band_hat)
            cand.valid = false;
            cand.reject_reason = 'no_shared_band_n0';
            CDResult = append_candidate(CDResult, cand);
            continue;
        end

        for m = 1:M
            if m > numel(y_tilde) || ~isfinite(y_tilde(m))
                cand.ap(m).valid = false; %#ok<AGROW>
                cand.ap(m).warning = 'invalid_ap_rss';
                continue;
            end

            ray = CD_project_prior_to_ray(q, APs.pos_xy(m, :), cfg);
            if ~ray.valid
                cand.ap(m).valid = false; %#ok<AGROW>
                cand.ap(m).warning = ray.reason;
                continue;
            end

            fit = build_shared_n0_ap_fit(y_tilde(m), P0, n_band_hat, ray, cfg);
            cand.ap(m).fit = fit; %#ok<AGROW>
            cand.ap(m).ray = ray;
            if ~fit.valid
                cand.ap(m).valid = false;
                cand.ap(m).warning = fit.warning;
                continue;
            end

            base_xy = make_base_points(q, grid_xy, ray, cfg);
            for k = 1:size(base_xy, 1)
                [val, w, patch_idx] = build_patch_prediction( ...
                    grid_xy, APs.pos_xy(m, :), base_xy(k, :), fit, cfg);
                if isempty(patch_idx)
                    continue;
                end
                num{b}(m, patch_idx) = num{b}(m, patch_idx) + cfg.calibration_weight * w(:)' .* val(:)';
                den{b}(m, patch_idx) = den{b}(m, patch_idx) + cfg.calibration_weight * w(:)';
            end
            cand.ap(m).valid = true;
            cand.ap(m).n_hat = fit.n_hat;
            cand.ap(m).PL0_hat_dB = fit.PL0_hat_dB;
            cand.ap(m).n0_source = ternary(m == ref_ap, 'nearest_ref_ap', 'band_shared');
        end

        cand.valid = any(arrayfun(@(x) isfield(x, 'valid') && x.valid, cand.ap));
        if cand.valid
            CDResult.n_events_valid = CDResult.n_events_valid + 1;
        else
            cand.reject_reason = 'no_valid_ap_fit';
        end
        CDResult = append_candidate(CDResult, cand);
    end

    for b = 1:B
        F_old = F_old_cache{b};
        F_new = F_old;
        touched = den{b} > cfg.blend_old_weight;
        pred = (cfg.blend_old_weight * F_old(touched) + num{b}(touched)) ./ max(den{b}(touched), eps);
        delta = pred - F_old(touched);
        delta = max(min(delta, cfg.max_delta_dB), -cfg.max_delta_dB);
        F_new(touched) = F_old(touched) + delta;

        PendingSpatialFP.band(b).F_dBm = F_new;
        PendingSpatialFP.band(b).RF_raw = F_new;
        PendingSpatialFP.band(b).F_lin = 10.^(F_new / 10);
        PendingSpatialFP.band(b).cd_uncertain_applied = any(touched(:));

        CDResult.band(b).n_entries_updated = nnz(touched);
        CDResult.band(b).mean_abs_delta_dB = mean_abs(F_new - F_old);
        CDResult.band(b).max_abs_delta_dB = max_abs(F_new - F_old);
        CDResult.band(b).status = ternary(any(touched(:)), 'ok', 'no_update');
    end

    PendingSpatialFP = CB_recompute_fingerprint_derivatives_from_dBm(PendingSpatialFP, Config);
    if CDResult.n_events_valid == 0
        CDResult.status = 'no_valid_uncertain_opportunistic_event';
    else
        CDResult.status = 'ok';
    end

    if cfg.verbose
        fprintf('[CD] uncertain opportunistic calibration: valid %d / %d events\n', ...
            CDResult.n_events_valid, CDResult.n_events_total);
        for b = 1:B
            fprintf('  band %d: %s, updated=%d, mean_abs_delta=%.3f dB, n0_median=%.3f\n', ...
                b, CDResult.band(b).status, CDResult.band(b).n_entries_updated, ...
                CDResult.band(b).mean_abs_delta_dB, CDResult.band(b).n0_hat_median);
        end
    end
end


function R = init_result(SpatialFP)
    R = struct();
    R.status = 'not_run';
    R.n_events_total = 0;
    R.n_events_valid = 0;
    R.candidates = struct([]);
    for b = 1:SpatialFP.B
        R.band(b).status = 'not_processed';
        R.band(b).n_entries_updated = 0;
        R.band(b).mean_abs_delta_dB = NaN;
        R.band(b).max_abs_delta_dB = NaN;
        R.band(b).n0_hat_values = [];
        R.band(b).n0_hat_median = NaN;
    end
end


function R = append_candidate(R, cand)
    if isempty(R.candidates)
        R.candidates = cand;
    else
        R.candidates(end+1) = cand; %#ok<AGROW>
    end
end


function [ok, reason] = accept_event(ev, SpatialFP, cfg)
    ok = false;
    reason = '';
    if ~isfield(ev, 'band_id') || ev.band_id < 1 || ev.band_id > SpatialFP.B
        reason = 'invalid_band';
        return;
    end
    if strcmpi(cfg.source_filter, 'opportunistic') && ...
            isfield(ev, 'source_type_name') && ~isempty(ev.source_type_name) && ...
            ~strcmpi(ev.source_type_name, 'opportunistic')
        reason = 'not_opportunistic';
        return;
    end
    if ~isfield(ev, 'location_prior') || ~isstruct(ev.location_prior) || ~isfield(ev.location_prior, 'type')
        reason = 'missing_location_prior';
        return;
    end
    if ~any(strcmpi(ev.location_prior.type, cfg.supported_prior_types))
        reason = 'unsupported_location_prior';
        return;
    end
    if ~isfield(ev, 'rss_dBm_by_ap') || isempty(ev.rss_dBm_by_ap)
        reason = 'missing_rss_by_ap';
        return;
    end
    if ~isfinite(resolve_power(ev, ev.band_id))
        reason = 'missing_exact_power';
        return;
    end
    ok = true;
end


function cand = rejected_candidate(ev, reason)
    cand = candidate_base(ev);
    cand.valid = false;
    cand.reject_reason = reason;
end


function cand = candidate_base(ev)
    cand = struct();
    cand.valid = false;
    cand.reject_reason = '';
    cand.source_uid = get_field_default(ev, 'source_uid', '');
    cand.band_id = get_field_default(ev, 'band_id', NaN);
    cand.source_type_name = get_field_default(ev, 'source_type_name', '');
    cand.frame_range = [get_field_default(ev, 'start_frame', NaN), get_field_default(ev, 'end_frame', NaN)];
    cand.P0_dBm = NaN;
    cand.P_e_dBm = NaN;
    cand.location_prior_type = '';
    cand.n0_mode = '';
    cand.n_band_hat = NaN;
    cand.n0_ref_ap = NaN;
    cand.n0_ref_fit = struct();
    cand.n0_ref_ray = struct();
    cand.ap = struct([]);
end


function [n0_by_band, n0_values] = estimate_shared_n0_by_band( ...
    calibration_events, SpatialFP, F_old_cache, APs, grid_xy, cfg)
    B = SpatialFP.B;
    n0_values = cell(1, B);
    n0_by_band = NaN(1, B);
    for b = 1:B
        n0_values{b} = [];
    end

    for e = 1:numel(calibration_events)
        ev = calibration_events(e);
        [ok_ev, ~] = accept_event(ev, SpatialFP, cfg);
        if ~ok_ev
            continue;
        end

        b = ev.band_id;
        P0 = get_offline_power(SpatialFP, b);
        P_e = resolve_power(ev, b);
        y_tilde = ev.rss_dBm_by_ap(:)' - (P_e - P0);
        q = CB_make_location_distribution_from_prior(ev.location_prior);

        [n_hat, ~, ~, ~, ~] = estimate_band_n0( ...
            q, F_old_cache{b}, y_tilde, P0, APs, grid_xy, cfg);
        if isfinite(n_hat)
            n0_values{b}(end+1) = n_hat; %#ok<AGROW>
        end
    end

    for b = 1:B
        if ~isempty(n0_values{b})
            n0_by_band(b) = median(n0_values{b});
        end
    end
end


function [n_hat, ref_ap, ref_fit, ref_ray, reason] = estimate_band_n0( ...
    q, F_old, y_tilde, P0, APs, grid_xy, cfg)
    n_hat = NaN;
    ref_ap = NaN;
    ref_fit = struct();
    ref_ray = struct();
    reason = '';

    center_xy = prior_center_xy(q);
    if any(~isfinite(center_xy))
        reason = 'invalid_prior_center_for_n0';
        return;
    end

    d_ap = sqrt(sum((APs.pos_xy - center_xy).^2, 2));
    [~, ref_ap] = min(d_ap);
    if ref_ap > numel(y_tilde) || ~isfinite(y_tilde(ref_ap))
        reason = 'nearest_ap_has_invalid_rss';
        return;
    end

    ref_ray = CD_project_prior_to_ray(q, APs.pos_xy(ref_ap, :), cfg);
    if ~ref_ray.valid
        reason = sprintf('nearest_ap_ray_invalid_%s', ref_ray.reason);
        return;
    end

    [r_old, y_old] = select_old_ray_samples( ...
        grid_xy, F_old(ref_ap, :)', APs.pos_xy(ref_ap, :), ref_ray, cfg);
    ref_fit = CD_fit_ray_lognormal_profile(r_old, y_old, P0, y_tilde(ref_ap), ref_ray, cfg);
    if ~ref_fit.valid
        reason = sprintf('nearest_ap_n0_fit_failed_%s', ref_fit.warning);
        return;
    end

    n_hat = ref_fit.n_hat;
end


function ref_ap = nearest_ap_for_prior(q, APs)
    center_xy = prior_center_xy(q);
    if any(~isfinite(center_xy))
        ref_ap = NaN;
        return;
    end
    d_ap = sqrt(sum((APs.pos_xy - center_xy).^2, 2));
    [~, ref_ap] = min(d_ap);
end


function fit = build_shared_n0_ap_fit(y_event_dBm, P0_dBm, n_hat, ray, cfg)
    fit = struct();
    fit.valid = false;
    fit.warning = 'not_run';
    fit.A_dBm = NaN;
    fit.B_dB_per_log10m = NaN;
    fit.n_hat = NaN;
    fit.PL0_hat_dB = NaN;
    fit.r_event_hat = NaN;
    fit.sigma_f_dB = NaN;
    fit.ell_m = NaN;
    fit.nmll = NaN;

    if ~isfinite(y_event_dBm) || ~isfinite(n_hat) || ~isfield(ray, 'mu_r') || ~isfinite(ray.mu_r)
        fit.warning = 'invalid_shared_n0_input';
        return;
    end

    r_anchor = max(ray.mu_r, cfg.eps_dist);
    B = -10 * n_hat;
    A = y_event_dBm - B * log10(r_anchor);

    fit.valid = true;
    fit.warning = '';
    fit.A_dBm = A;
    fit.B_dB_per_log10m = B;
    fit.n_hat = n_hat;
    fit.PL0_hat_dB = P0_dBm - A;
    fit.r_event_hat = r_anchor;
    fit.sigma_f_dB = NaN;
    fit.ell_m = NaN;
    fit.nmll = NaN;
end


function xy = prior_center_xy(q)
    switch lower(q.type)
        case 'rect'
            bbox = q.bbox;
            xy = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
        case 'gaussian'
            xy = q.mu(:)';
        case 'point'
            xy = q.x(:)';
        otherwise
            xy = [NaN NaN];
    end
end


function [r_old, y_old] = select_old_ray_samples(grid_xy, y_map, ap_xy, ray, cfg)
    rel = grid_xy - ap_xy;
    r = rel * ray.dir(:);
    perp = sqrt(max(sum(rel.^2, 2) - r.^2, 0));
    mask = r > cfg.eps_dist & perp <= cfg.ray_width_m & isfinite(y_map);
    if sum(mask) < 8
        d_center = sqrt(sum((grid_xy - ray.center_xy).^2, 2));
        [~, ord] = sort(d_center, 'ascend');
        ord = ord(1:min(numel(ord), cfg.max_ray_points));
        r_old = sqrt(sum((grid_xy(ord, :) - ap_xy).^2, 2));
        y_old = y_map(ord);
        return;
    end
    r_old = r(mask);
    y_old = y_map(mask);
    [r_old, ord] = sort(r_old, 'ascend');
    y_old = y_old(ord);
end


function base_xy = make_base_points(q, grid_xy, ray, cfg)
    switch lower(q.type)
        case 'rect'
            bbox = q.bbox;
            mask = grid_xy(:,1) >= bbox(1) & grid_xy(:,1) <= bbox(2) & ...
                   grid_xy(:,2) >= bbox(3) & grid_xy(:,2) <= bbox(4);
            pts = grid_xy(mask, :);
        case 'gaussian'
            S = 0.5 * (q.Sigma + q.Sigma') + 1e-6 * eye(2);
            X = grid_xy - q.mu;
            maha = sum((X / S) .* X, 2);
            pts = grid_xy(maha <= cfg.gaussian_conf_chi2, :);
        otherwise
            pts = [];
    end
    if isempty(pts)
        pts = ray.center_xy;
    end
    d = sqrt(sum((pts - ray.center_xy).^2, 2));
    [~, ord] = sort(d, 'ascend');
    ord = ord(1:min(numel(ord), cfg.max_base_points));
    base_xy = pts(ord, :);
end


function [val, w, idx] = build_patch_prediction(grid_xy, ap_xy, base_xy, fit, cfg)
    d_patch = sqrt(sum((grid_xy - base_xy).^2, 2));
    idx = find(d_patch <= cfg.patch_radius_m);
    if isempty(idx)
        val = [];
        w = [];
        return;
    end
    d_base = max(norm(base_xy - ap_xy), cfg.eps_dist);
    y_base = fit.A_dBm + fit.B_dB_per_log10m * log10(d_base);
    d_grid = sqrt(sum((grid_xy(idx, :) - ap_xy).^2, 2));
    d_grid = max(d_grid, cfg.eps_dist);
    val = y_base - 10 * fit.n_hat * log10(d_grid / d_base);

    taper_width = min(max(cfg.patch_taper_width_m, 0), cfg.patch_radius_m);
    w = ones(numel(idx), 1);
    if taper_width > 0
        edge0 = cfg.patch_radius_m - taper_width;
        edge = d_patch(idx) > edge0;
        z = (d_patch(idx(edge)) - edge0) / taper_width;
        w(edge) = 0.5 * (1 + cos(pi * z));
    end
end


function P = resolve_power(ev, band_id)
    P = NaN;
    if isfield(ev, 'P_e_dBm') && isfinite(ev.P_e_dBm)
        P = ev.P_e_dBm;
        return;
    end
    if ~isfield(ev, 'power_prior') || ~isstruct(ev.power_prior) || ~isfield(ev.power_prior, 'type')
        return;
    end
    pp = ev.power_prior;
    switch lower(pp.type)
        case 'exact'
            if isscalar(pp.value)
                P = pp.value;
            elseif band_id <= numel(pp.value)
                P = pp.value(band_id);
            end
        case 'exact_by_band'
            if band_id <= numel(pp.value)
                P = pp.value(band_id);
            end
    end
end


function P0 = get_offline_power(SpatialFP, b)
    if isfield(SpatialFP.band(b), 'reference_power_dBm_used')
        P0 = SpatialFP.band(b).reference_power_dBm_used;
    elseif isfield(SpatialFP, 'reference_power_dBm_used')
        P = SpatialFP.reference_power_dBm_used;
        P0 = P(min(b, numel(P)));
    elseif isfield(SpatialFP, 'ref_power_dBm')
        P = SpatialFP.ref_power_dBm;
        P0 = P(min(b, numel(P)));
    else
        error('[CD] missing offline reference power for band %d', b);
    end
end


function F = get_band_dbm(band_fp)
    if isfield(band_fp, 'F_dBm')
        F = band_fp.F_dBm;
    elseif isfield(band_fp, 'RF_raw')
        F = band_fp.RF_raw;
    else
        error('[CD] band fingerprint missing F_dBm/RF_raw');
    end
end


function v = mean_abs(X)
    x = X(:);
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = mean(abs(x));
    end
end


function v = max_abs(X)
    x = X(:);
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = max(abs(x));
    end
end


function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end


function v = get_field_default(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end
