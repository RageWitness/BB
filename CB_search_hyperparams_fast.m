function [Config_best, SearchResult] = CB_search_hyperparams_fast( ...
    SpatialFP_current, calibration_events, validation_events, FrameStates, Config)
% CB_SEARCH_HYPERPARAMS_FAST  Fast ell/trust-ratio search for CB calibration.
%
% Search uses a cheap Gaussian-kernel RSS correction map. The selected
% parameters are intended for one final full CB-GPR update.

    if nargin < 5
        Config = struct();
    end

    cfg = CB_calib_ap_profile_defaults(Config);
    Config_best = Config;
    Config_best.calib_ap_profile = cfg;

    SearchResult = struct();
    SearchResult.status = 'not_run';
    SearchResult.backend = cfg.fast_search_backend;
    SearchResult.table = struct([]);
    SearchResult.best_ell_m = cfg.ell_x_m;
    SearchResult.best_trust_ratio = cfg.sigma0_dB / max(cfg.sigma_meas_dB, eps);
    SearchResult.best_rmse_m = NaN;

    if isempty(validation_events)
        SearchResult.status = 'no_validation_events';
        return;
    end
    if isempty(calibration_events)
        SearchResult.status = 'no_calibration_events';
        return;
    end

    ell_grid = cfg.fast_ell_grid_m(:)';
    trust_grid = cfg.fast_trust_ratio_grid(:)';
    [best_ell, best_ratio, best_rmse, table1] = run_grid( ...
        SpatialFP_current, calibration_events, validation_events, FrameStates, Config, ell_grid, trust_grid);

    table_all = table1;
    if cfg.fast_refine_enable && isfinite(best_rmse)
        ell_ref = unique(max(1, best_ell * cfg.fast_refine_scale(:)'));
        ratio_ref = unique(max(0.25, best_ratio * cfg.fast_refine_scale(:)'));
        [best_ell2, best_ratio2, best_rmse2, table2] = run_grid( ...
            SpatialFP_current, calibration_events, validation_events, FrameStates, Config, ell_ref, ratio_ref);
        table_all = [table_all, table2]; %#ok<AGROW>
        if best_rmse2 <= best_rmse
            best_ell = best_ell2;
            best_ratio = best_ratio2;
            best_rmse = best_rmse2;
        end
    end

    cfg.ell_x_m = best_ell;
    cfg.ell_y_m = best_ell;
    cfg.sigma0_dB = best_ratio * cfg.sigma_meas_dB;
    Config_best.calib_ap_profile = cfg;

    SearchResult.status = 'ok';
    SearchResult.table = table_all;
    SearchResult.best_ell_m = best_ell;
    SearchResult.best_trust_ratio = best_ratio;
    SearchResult.best_rmse_m = best_rmse;

    fprintf('[CB-fast-search] best ell=%.2f m, trust_ratio=%.2f, validation RMSE=%.3f m\n', ...
        best_ell, best_ratio, best_rmse);
end


function [best_ell, best_ratio, best_rmse, rows] = run_grid( ...
    SpatialFP_current, calibration_events, validation_events, FrameStates, Config, ell_grid, trust_grid)

    cfg0 = CB_calib_ap_profile_defaults(Config);
    rows = struct([]);
    best_rmse = inf;
    best_ell = ell_grid(1);
    best_ratio = trust_grid(1);
    n = 0;

    for ie = 1:numel(ell_grid)
        for ir = 1:numel(trust_grid)
            ell = ell_grid(ie);
            ratio = trust_grid(ir);
            cfg = cfg0;
            cfg.ell_x_m = ell;
            cfg.ell_y_m = ell;
            cfg.sigma0_dB = ratio * cfg.sigma_meas_dB;

            Config_i = Config;
            Config_i.calib_ap_profile = cfg;

            tic;
            SpatialFP_fast = build_fast_pending_map(SpatialFP_current, calibration_events, cfg, Config_i);
            rmse = compute_validation_rmse(validation_events, SpatialFP_fast, FrameStates, Config_i);
            elapsed = toc;

            n = n + 1;
            rows(n).ell_m = ell; %#ok<AGROW>
            rows(n).trust_ratio = ratio;
            rows(n).sigma0_dB = cfg.sigma0_dB;
            rows(n).sigma_meas_dB = cfg.sigma_meas_dB;
            rows(n).rmse_m = rmse;
            rows(n).elapsed_s = elapsed;

            fprintf('[CB-fast-search] ell=%6.2f ratio=%5.2f rmse=%8.3f m (%.2fs)\n', ...
                ell, ratio, rmse, elapsed);

            if rmse < best_rmse
                best_rmse = rmse;
                best_ell = ell;
                best_ratio = ratio;
            end
        end
    end
end


function SpatialFP_fast = build_fast_pending_map(SpatialFP_current, events, cfg, Config)
    SpatialFP_fast = SpatialFP_current;
    grid_xy = SpatialFP_current.grid_xy;
    B = SpatialFP_current.B;
    M = SpatialFP_current.M;

    for b = 1:B
        F_old = SpatialFP_current.band(b).F_dBm;
        F_new = F_old;
        idx_b = find([events.band_id] == b);
        if isempty(idx_b)
            continue;
        end

        for m = 1:M
            num = zeros(1, SpatialFP_current.G);
            den = zeros(1, SpatialFP_current.G);
            for kk = 1:numel(idx_b)
                ev = events(idx_b(kk));
                [ok, y_tilde] = normalized_event_rss(ev, b, m, SpatialFP_current);
                if ~ok || ~isfinite(y_tilde)
                    continue;
                end
                q = CB_make_location_distribution_from_prior(ev.location_prior);
                x0 = representative_xy(q);
                y0 = profile_value_at_q(grid_xy, F_old(m, :)', q);
                residual = y_tilde - y0;
                nu = event_noise_var(ev, grid_xy, F_old(m, :)', q, b, m, cfg);
                w = 1 / max(nu, eps);

                d2 = sum((grid_xy - x0).^2, 2)';
                k = exp(-0.5 * d2 / max(cfg.ell_x_m^2, eps));
                num = num + w * residual * k;
                den = den + w * k;
            end

            old_w = 1 / max(cfg.sigma0_dB^2, eps);
            delta = num ./ (old_w + den);
            F_new(m, :) = F_old(m, :) + delta;
        end

        SpatialFP_fast.band(b).F_dBm = F_new;
    end

    SpatialFP_fast = CB_recompute_fingerprint_derivatives_from_dBm(SpatialFP_fast, Config);
    SpatialFP_fast.is_pending_cb_fast_map = true;
end


function rmse = compute_validation_rmse(validation_events, SpatialFP_fast, FrameStates, Config)
    ev = validation_events;
    for i = 1:numel(ev)
        ev(i).route_action = 'localize_only';
        if ~isfield(ev(i), 'type_hat') || isempty(ev(i).type_hat)
            ev(i).type_hat = 'cb_validation';
        end
    end
    LocResults = run_m4_wknn_localization(ev, SpatialFP_fast, FrameStates, Config);
    if isempty(LocResults) || ~isfield(LocResults, 'loc_error')
        rmse = inf;
        return;
    end
    e = [LocResults.loc_error];
    e = e(isfinite(e));
    if isempty(e)
        rmse = inf;
    else
        rmse = sqrt(mean(e.^2));
    end
end


function [ok, y_tilde] = normalized_event_rss(ev, band_id, ap_id, SpatialFP)
    ok = false;
    y_tilde = NaN;
    y = event_rss_by_ap(ev);
    if numel(y) < ap_id || ~isfinite(y(ap_id))
        return;
    end
    [ok_p, Pe] = event_power_dBm(ev, band_id);
    if ~ok_p
        return;
    end
    P0 = offline_power_dBm(SpatialFP, band_id);
    y_tilde = y(ap_id) - (Pe - P0);
    ok = true;
end


function y = event_rss_by_ap(ev)
    if isfield(ev, 'rss_dBm_by_ap') && ~isempty(ev.rss_dBm_by_ap)
        y = ev.rss_dBm_by_ap(:)';
    elseif isfield(ev, 'obs_segment_dBm') && ~isempty(ev.obs_segment_dBm)
        X = ev.obs_segment_dBm;
        if isvector(X)
            y = X(:)';
        else
            y = nan(1, size(X, 1));
            for m = 1:size(X, 1)
                v = X(m, :);
                v = v(isfinite(v));
                if ~isempty(v), y(m) = mean(v); end
            end
        end
    else
        y = [];
    end
end


function [ok, P] = event_power_dBm(ev, band_id)
    ok = false;
    P = NaN;
    if isfield(ev, 'P_e_dBm') && isfinite(ev.P_e_dBm)
        P = ev.P_e_dBm;
        ok = true;
        return;
    end
    if ~isfield(ev, 'power_prior') || ~isstruct(ev.power_prior) || ~isfield(ev.power_prior, 'type')
        return;
    end
    pp = ev.power_prior;
    switch pp.type
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
    ok = isfinite(P);
end


function P0 = offline_power_dBm(SpatialFP, band_id)
    if isfield(SpatialFP.band(band_id), 'reference_power_dBm_used')
        P0 = SpatialFP.band(band_id).reference_power_dBm_used;
    elseif isfield(SpatialFP, 'reference_power_dBm_used')
        P = SpatialFP.reference_power_dBm_used;
        P0 = P(min(band_id, numel(P)));
    elseif isfield(SpatialFP, 'ref_power_dBm')
        P = SpatialFP.ref_power_dBm;
        P0 = P(min(band_id, numel(P)));
    else
        error('[CB-fast-search] missing offline reference power');
    end
end


function x = representative_xy(q)
    switch q.type
        case 'point'
            x = q.x;
        case 'rect'
            x = [0.5 * (q.bbox(1) + q.bbox(2)), 0.5 * (q.bbox(3) + q.bbox(4))];
        case 'gaussian'
            x = q.mu;
        case 'trajectory'
            [pts, w] = CB_trajectory_quadrature(q, 15);
            x = sum(pts .* w, 1);
        otherwise
            x = [0, 0];
    end
end


function y = profile_value_at_q(grid_xy, f, q)
    switch q.type
        case 'point'
            y = nearest_value(grid_xy, f, q.x);
        case 'rect'
            mask = grid_xy(:,1) >= q.bbox(1) & grid_xy(:,1) <= q.bbox(2) & ...
                   grid_xy(:,2) >= q.bbox(3) & grid_xy(:,2) <= q.bbox(4);
            vals = f(mask);
            vals = vals(isfinite(vals));
            if isempty(vals)
                y = nearest_value(grid_xy, f, representative_xy(q));
            else
                y = mean(vals);
            end
        case 'gaussian'
            d2 = sum((grid_xy - q.mu).^2, 2);
            sig2 = max(mean(diag(q.Sigma)), 1);
            w = exp(-0.5 * d2 / sig2);
            y = sum(w .* f) / max(sum(w), eps);
        case 'trajectory'
            [pts, w] = CB_trajectory_quadrature(q, 15);
            vals = zeros(numel(w), 1);
            for i = 1:numel(w)
                vals(i) = nearest_value(grid_xy, f, pts(i, :));
            end
            y = sum(vals .* w);
        otherwise
            y = nearest_value(grid_xy, f, representative_xy(q));
    end
end


function y = nearest_value(grid_xy, f, x)
    d2 = sum((grid_xy - x).^2, 2);
    [~, idx] = min(d2);
    y = f(idx);
end


function nu = event_noise_var(ev, grid_xy, f, q, b, m, cfg)
    sigma = param_value(cfg.sigma_meas_dB, b, m);
    nu = sigma^2;
    if strcmp(q.type, 'rect')
        [v_region, ~] = CB_estimate_region_rss_variance(grid_xy, f, q, cfg);
        nu = nu + cfg.gamma_region_var * v_region;
    elseif strcmp(q.type, 'gaussian') && cfg.gamma_gaussian_var > 0
        y0 = profile_value_at_q(grid_xy, f, q);
        d2 = sum((grid_xy - q.mu).^2, 2);
        sig2 = max(mean(diag(q.Sigma)), 1);
        w = exp(-0.5 * d2 / sig2);
        nu = nu + cfg.gamma_gaussian_var * sum(w .* (f - y0).^2) / max(sum(w), eps);
    end
    if isfield(ev, 'location_prior') && isfield(ev.location_prior, 'type') && strcmp(ev.location_prior.type, 'region')
        nu = max(nu, sigma^2);
    end
end


function v = param_value(x, b, m)
    if isscalar(x)
        v = x;
    elseif ismatrix(x) && size(x, 1) >= b && size(x, 2) >= m
        v = x(b, m);
    elseif isvector(x) && numel(x) >= m
        v = x(m);
    elseif isvector(x) && numel(x) >= b
        v = x(b);
    else
        v = x(1);
    end
end
