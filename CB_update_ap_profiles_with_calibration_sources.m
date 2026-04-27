function [APProfile_pending, SpatialFP_pending, diagnostics] = CB_update_ap_profiles_with_calibration_sources( ...
    APProfile_current, SpatialFP_current, calibration_events, Config)
% CB_UPDATE_AP_PROFILES_WITH_CALIBRATION_SOURCES  Update AP profiles and RSS dBm library.

    if nargin < 4
        Config = struct();
    end
    cfg = CB_calib_ap_profile_defaults(Config);

    APProfile_pending = APProfile_current;
    SpatialFP_pending = SpatialFP_current;
    diagnostics = init_diagnostics(SpatialFP_current);
    diagnostics.config = cfg;

    if isempty(calibration_events)
        diagnostics.status = 'no_calibration_events';
        return;
    end

    B = SpatialFP_current.B;
    M = SpatialFP_current.M;
    G = SpatialFP_current.G;
    grid_xy = SpatialFP_current.grid_xy;

    event_bank = prepare_events(calibration_events, SpatialFP_current, cfg);
    diagnostics.n_events_total = numel(calibration_events);
    diagnostics.n_events_valid = numel(event_bank);
    if isempty(event_bank)
        diagnostics.status = 'no_valid_power_known_events';
        return;
    end

    for b = 1:B
        F_old = get_band_dbm(SpatialFP_current.band(b));
        F_new = F_old;
        ev_idx_b = find([event_bank.band_id] == b);
        if isempty(ev_idx_b)
            diagnostics.band(b).status = 'no_events_for_band';
            continue;
        end

        for m = 1:M
            ap = APProfile_current.band(b).ap(m);
            kp = ap.kernel_params;
            train = build_training_set_for_ap(ap, event_bank(ev_idx_b), grid_xy, F_old(m, :)', b, m, cfg);

            [mu_pred, var_pred, gp_info] = predict_profile_grid(train, grid_xy, ap.mean_model, kp, cfg);
            F_new(m, :) = mu_pred(:)';

            APProfile_pending.band(b).ap(m).train_inputs = train.q;
            APProfile_pending.band(b).ap(m).train_outputs = train.Y;
            APProfile_pending.band(b).ap(m).train_noise = train.nu;
            APProfile_pending.band(b).ap(m).posterior_mean_dBm = mu_pred(:)';
            APProfile_pending.band(b).ap(m).posterior_var_dB2 = var_pred(:)';
            APProfile_pending.band(b).ap(m).update_info = gp_info;

            delta = F_new(m, :) - F_old(m, :);
            diagnostics.band(b).ap(m).max_abs_delta = max(abs(delta(isfinite(delta))));
            diagnostics.band(b).ap(m).mean_abs_delta = mean_abs(delta);
            diagnostics.band(b).ap(m).rmse_delta = sqrt(mean(delta(isfinite(delta)).^2));
            diagnostics.band(b).ap(m).n_train = numel(train.Y);
            diagnostics.band(b).ap(m).n_calibration_train = train.n_calibration;
            diagnostics.band(b).ap(m).backend = gp_info.backend;
        end

        SpatialFP_pending.band(b).F_dBm = F_new;
        diagnostics.band(b).status = 'ok';
        diagnostics.band(b).mean_abs_delta = mean_abs(F_new - F_old);
        diagnostics.band(b).max_abs_delta = max_abs(F_new - F_old);
    end

    SpatialFP_pending = CB_recompute_fingerprint_derivatives_from_dBm(SpatialFP_pending, Config);
    SpatialFP_pending.is_pending_cb_map = true;
    SpatialFP_pending.cb_method = 'ap_profile_uncertain_input_gp';
    APProfile_pending.update_round = get_field_default(APProfile_current, 'update_round', 0) + 1;
    APProfile_pending.last_update_diagnostics = diagnostics;
    diagnostics.status = 'ok';

    if cfg.verbose
        fprintf('[CB] AP-profile update: valid events %d / %d\n', ...
            diagnostics.n_events_valid, diagnostics.n_events_total);
        for b = 1:B
            fprintf('  band %d: %s, mean_abs_delta=%.3f dB\n', ...
                b, diagnostics.band(b).status, diagnostics.band(b).mean_abs_delta);
        end
    end
end


function diagnostics = init_diagnostics(SpatialFP)
    diagnostics = struct();
    diagnostics.status = 'not_run';
    diagnostics.n_events_total = 0;
    diagnostics.n_events_valid = 0;
    for b = 1:SpatialFP.B
        diagnostics.band(b).status = 'not_processed';
        diagnostics.band(b).mean_abs_delta = NaN;
        diagnostics.band(b).max_abs_delta = NaN;
        for m = 1:SpatialFP.M
            diagnostics.band(b).ap(m).max_abs_delta = NaN;
            diagnostics.band(b).ap(m).mean_abs_delta = NaN;
            diagnostics.band(b).ap(m).rmse_delta = NaN;
            diagnostics.band(b).ap(m).n_train = 0;
            diagnostics.band(b).ap(m).n_calibration_train = 0;
            diagnostics.band(b).ap(m).backend = '';
        end
    end
end


function events = prepare_events(calibration_events, SpatialFP, cfg)
    events = struct([]);
    n = 0;
    for i = 1:numel(calibration_events)
        ev0 = calibration_events(i);
        if ~isfield(ev0, 'band_id') || ev0.band_id < 1 || ev0.band_id > SpatialFP.B
            continue;
        end
        [ok_power, P_e] = resolve_power(ev0, ev0.band_id);
        if ~ok_power
            if cfg.reject_power_none
                fprintf('[CB] reject event %d: current AP_profile calibration requires known power\n', i);
            end
            continue;
        end
        q = CB_make_location_distribution_from_prior(ev0.location_prior);
        P0 = get_offline_power(SpatialFP, ev0.band_id);
        y = get_event_rss_by_ap(ev0);
        if isempty(y)
            continue;
        end
        y_tilde = y(:)' - (P_e - P0);

        ev = struct();
        ev.source_uid = get_field_default(ev0, 'source_uid', sprintf('cb_event_%d', i));
        ev.band_id = ev0.band_id;
        ev.label = get_field_default(ev0, 'label', []);
        ev.source_type_name = get_field_default(ev0, 'source_type_name', '');
        ev.start_frame = get_field_default(ev0, 'start_frame', NaN);
        ev.end_frame = get_field_default(ev0, 'end_frame', NaN);
        ev.time_range = get_field_default(ev0, 'time_range', [ev.start_frame, ev.end_frame]);
        ev.location_prior = ev0.location_prior;
        ev.power_prior = ev0.power_prior;
        ev.metadata = get_field_default(ev0, 'metadata', struct());
        ev.template_key = get_field_default(ev0, 'template_key', '');
        ev.q = q;
        ev.P_e_dBm = P_e;
        ev.P0_dBm = P0;
        ev.rss_dBm_by_ap = y(:)';
        ev.rss_tilde_dBm_by_ap = y_tilde;
        ev.power_normalization_delta_dB = P_e - P0;
        n = n + 1;
        if n == 1
            events = ev;
        else
            events(n) = ev; %#ok<AGROW>
        end
    end
end


function train = build_training_set_for_ap(ap, events_b, grid_xy, F_old_m, b, m, cfg)
    idx = select_offline_subset(grid_xy, events_b, cfg.max_offline_train_points);
    n0 = numel(idx);
    ncal = numel(events_b);

    q = cell(n0 + ncal, 1);
    Y = zeros(n0 + ncal, 1);
    nu = zeros(n0 + ncal, 1);
    source = cell(n0 + ncal, 1);

    for i = 1:n0
        q{i} = struct('type', 'point', 'x', grid_xy(idx(i), :));
        Y(i) = F_old_m(idx(i));
        nu(i) = get_scalar_or_vec(ap.offline_noise_var, idx(i));
        source{i} = 'offline';
    end

    for e = 1:ncal
        row = n0 + e;
        q{row} = events_b(e).q;
        Y(row) = events_b(e).rss_tilde_dBm_by_ap(m);
        nu(row) = measurement_noise_for_event(events_b(e), grid_xy, F_old_m, m, cfg);
        source{row} = 'calibration';
    end

    valid = isfinite(Y) & isfinite(nu) & nu > 0;
    train = struct();
    train.q = q(valid);
    train.Y = Y(valid);
    train.nu = nu(valid);
    train.source = source(valid);
    train.n_offline = sum(strcmp(train.source, 'offline'));
    train.n_calibration = sum(strcmp(train.source, 'calibration'));
    train.band_id = b;
    train.ap_id = m;
end


function idx = select_offline_subset(grid_xy, events_b, max_n)
    G = size(grid_xy, 1);
    max_n = min(max_n, G);
    if G <= max_n
        idx = (1:G)';
        return;
    end

    seed_pts = [];
    for e = 1:numel(events_b)
        seed_pts = [seed_pts; representative_xy(events_b(e).q)]; %#ok<AGROW>
    end
    if isempty(seed_pts)
        ord = round(linspace(1, G, max_n));
        idx = unique(ord(:), 'stable');
        return;
    end

    dmin = inf(G, 1);
    for s = 1:size(seed_pts, 1)
        d2 = sum((grid_xy - seed_pts(s, :)).^2, 2);
        dmin = min(dmin, d2);
    end
    [~, near_ord] = sort(dmin, 'ascend');
    n_near = min(ceil(0.6 * max_n), G);
    near_idx = near_ord(1:n_near);
    stride_idx = round(linspace(1, G, max_n - n_near));
    idx = unique([near_idx(:); stride_idx(:)], 'stable');
    if numel(idx) > max_n
        idx = idx(1:max_n);
    end
end


function xy = representative_xy(q)
    switch lower(q.type)
        case 'point'
            xy = q.x;
        case 'rect'
            xy = [0.5 * (q.bbox(1) + q.bbox(2)), 0.5 * (q.bbox(3) + q.bbox(4))];
        case 'gaussian'
            xy = q.mu;
        case 'trajectory'
            [pts, w] = CB_trajectory_quadrature(q, 15);
            xy = sum(pts .* w, 1);
        otherwise
            xy = [0, 0];
    end
end


function nu = measurement_noise_for_event(ev, grid_xy, F_old_m, ap_id, cfg)
    sigma_meas = get_param_value(cfg.sigma_meas_dB, ev.band_id, ap_id);
    nu = sigma_meas^2;
    if strcmp(ev.q.type, 'rect')
        [v_region, info] = CB_estimate_region_rss_variance(grid_xy, F_old_m, ev.q, cfg);
        nu = nu + cfg.gamma_region_var * v_region;
        ev.region_variance_info = info; %#ok<NASGU>
    elseif strcmp(ev.q.type, 'gaussian') && cfg.gamma_gaussian_var > 0
        v_gauss = estimate_gaussian_var(grid_xy, F_old_m, ev.q, cfg);
        nu = nu + cfg.gamma_gaussian_var * v_gauss;
    end
end


function v = estimate_gaussian_var(grid_xy, F_old_m, q, cfg)
    d2 = sum((grid_xy - q.mu).^2, 2);
    k = min(numel(d2), max(10, cfg.gaussian_quad_n^2));
    [~, ord] = sort(d2, 'ascend');
    vals = F_old_m(ord(1:k));
    vals = vals(isfinite(vals));
    if numel(vals) <= 1
        v = 0;
    else
        v = var(vals, 0);
    end
end


function [mu_pred, var_pred, info] = predict_profile_grid(train, grid_xy, mean_model, kp, cfg)
    N = numel(train.Y);
    if N == 0
        mu_pred = CB_eval_ap_profile_mean(mean_model, grid_xy);
        var_pred = kp.sigma_f_dB^2 * ones(size(mu_pred));
        info = struct('backend', 'mean_only', 'n_train', 0);
        return;
    end

    if N > cfg.max_train_points_exact_gp
        error('[CB] training set has %d points, above max_train_points_exact_gp=%d after scalable subset', ...
            N, cfg.max_train_points_exact_gp);
    end

    K = zeros(N, N);
    for i = 1:N
        for j = i:N
            kij = CB_compute_uncertain_kernel(train.q{i}, train.q{j}, kp);
            K(i, j) = kij;
            K(j, i) = kij;
        end
    end
    K = 0.5 * (K + K') + diag(train.nu(:));
    jitter = max(cfg.jitter_abs, cfg.jitter_rel * mean(max(diag(K), eps)));
    [L, jitter_used] = stable_chol(K, jitter);

    mu_train = eval_mean_for_q_set(mean_model, train.q, cfg);
    alpha = L' \ (L \ (train.Y(:) - mu_train(:)));

    G = size(grid_xy, 1);
    mu_pred = zeros(G, 1);
    var_pred = zeros(G, 1);
    block = max(1, cfg.prediction_block_size);
    for st = 1:block:G
        ed = min(G, st + block - 1);
        xb = grid_xy(st:ed, :);
        nb = size(xb, 1);
        Kstar = zeros(nb, N);
        for ii = 1:nb
            qx = struct('type', 'point', 'x', xb(ii, :));
            for j = 1:N
                Kstar(ii, j) = CB_compute_uncertain_kernel(qx, train.q{j}, kp);
            end
        end
        mu0 = CB_eval_ap_profile_mean(mean_model, xb);
        mu_pred(st:ed) = mu0 + Kstar * alpha;
        V = L \ Kstar';
        var_pred(st:ed) = max(kp.sigma_f_dB^2 - sum(V.^2, 1)', 0);
    end

    info = struct();
    info.backend = 'scalable_subset_gp';
    info.n_train = N;
    info.jitter_used = jitter_used;
end


function mu = eval_mean_for_q_set(mean_model, q, cfg)
    mu = zeros(numel(q), 1);
    for i = 1:numel(q)
        mu(i) = eval_mean_for_q(mean_model, q{i}, cfg);
    end
end


function mu = eval_mean_for_q(mean_model, q, cfg)
    switch lower(q.type)
        case 'point'
            mu = CB_eval_ap_profile_mean(mean_model, q.x);
        case 'rect'
            n = max(3, cfg.rect_quad_n);
            xs = linspace(q.bbox(1), q.bbox(2), n);
            ys = linspace(q.bbox(3), q.bbox(4), n);
            [XX, YY] = meshgrid(xs, ys);
            mu = mean(CB_eval_ap_profile_mean(mean_model, [XX(:), YY(:)]));
        case 'gaussian'
            pts = gaussian_sigma_points(q.mu, q.Sigma);
            mu = mean(CB_eval_ap_profile_mean(mean_model, pts));
        case 'trajectory'
            [pts, w] = CB_trajectory_quadrature(q, cfg.trajectory_quad_n);
            mu = sum(CB_eval_ap_profile_mean(mean_model, pts) .* w);
        otherwise
            mu = 0;
    end
end


function pts = gaussian_sigma_points(mu, Sigma)
    mu = mu(:)';
    Sigma = 0.5 * (Sigma + Sigma');
    [V, D] = eig(Sigma);
    s = sqrt(max(diag(D), 0))';
    pts = [mu; mu + s(1) * V(:,1)'; mu - s(1) * V(:,1)'; ...
               mu + s(2) * V(:,2)'; mu - s(2) * V(:,2)'];
end


function [L, jitter_used] = stable_chol(K, jitter)
    jitter_used = jitter;
    N = size(K, 1);
    for a = 1:8
        [L, p] = chol(K + jitter_used * eye(N), 'lower');
        if p == 0
            return;
        end
        jitter_used = jitter_used * 10;
    end
    error('[CB] Cholesky failed for AP-profile covariance');
end


function y = get_event_rss_by_ap(ev)
    if isfield(ev, 'rss_dBm_by_ap') && ~isempty(ev.rss_dBm_by_ap)
        y = ev.rss_dBm_by_ap(:)';
    elseif isfield(ev, 'obs_segment_dBm') && ~isempty(ev.obs_segment_dBm)
        X = ev.obs_segment_dBm;
        if isvector(X)
            y = X(:)';
        else
            y = finite_mean_rows(X);
        end
    else
        y = [];
    end
end


function y = finite_mean_rows(X)
    y = nan(1, size(X, 1));
    for i = 1:size(X, 1)
        v = X(i, :);
        v = v(isfinite(v));
        if ~isempty(v)
            y(i) = mean(v);
        end
    end
end


function [ok, P_e] = resolve_power(ev, band_id)
    ok = false;
    P_e = NaN;
    if isfield(ev, 'P_e_dBm') && isfinite(ev.P_e_dBm)
        P_e = ev.P_e_dBm;
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
                P_e = pp.value;
            elseif band_id <= numel(pp.value)
                P_e = pp.value(band_id);
            end
        case 'exact_by_band'
            if band_id <= numel(pp.value)
                P_e = pp.value(band_id);
            end
    end
    ok = isfinite(P_e);
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
        error('[CB] missing offline reference power for band %d', b);
    end
end


function F = get_band_dbm(band_fp)
    if isfield(band_fp, 'F_dBm')
        F = band_fp.F_dBm;
    elseif isfield(band_fp, 'RF_raw')
        F = band_fp.RF_raw;
    else
        error('[CB] band fingerprint missing F_dBm/RF_raw');
    end
end


function v = get_scalar_or_vec(x, idx)
    if isscalar(x)
        v = x;
    else
        v = x(min(idx, numel(x)));
    end
end


function v = get_param_value(x, b, m)
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


function v = get_field_default(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end
