function [APProfile, SpatialFP_init] = CB_build_ap_profiles_from_offline_library(SpatialFP, grid_xy, Config)
% CB_BUILD_AP_PROFILES_FROM_OFFLINE_LIBRARY  Build AP-profile model metadata.

    if nargin < 2 || isempty(grid_xy)
        grid_xy = SpatialFP.grid_xy;
    end
    if nargin < 3
        Config = struct();
    end

    cfg = CB_calib_ap_profile_defaults(Config);
    B = SpatialFP.B;
    M = SpatialFP.M;
    G = SpatialFP.G;

    APProfile = struct();
    APProfile.method = 'CB_ap_profile_uncertain_input_gp';
    APProfile.grid_xy = grid_xy;
    APProfile.B = B;
    APProfile.M = M;
    APProfile.G = G;
    APProfile.config = cfg;
    APProfile.created_from = 'offline_spatialfp';

    for b = 1:B
        F0 = get_band_dbm(SpatialFP.band(b));
        APProfile.band(b).status = 'ok';
        APProfile.band(b).band_id = b;
        APProfile.band(b).reference_power_dBm_used = get_offline_power(SpatialFP, b);

        for m = 1:M
            y0 = F0(m, :)';
            mean_model = fit_log_distance_mean(grid_xy, y0, m, SpatialFP, Config, b);
            kp = make_kernel_params(cfg, b, m);
            [kp, offline_fit] = calibrate_offline_profile_hyperparams( ...
                grid_xy, y0, mean_model, kp, cfg, b, m);

            ap = struct();
            ap.kernel_params = kp;
            ap.mean_model = mean_model;
            ap.offline_xy = grid_xy;
            ap.offline_y_dBm = y0;
            ap.offline_noise_var = offline_fit.sigma0_dB^2 * ones(G, 1);
            ap.offline_fit = offline_fit;
            ap.posterior_mean_dBm = y0(:)';
            ap.posterior_var_dB2 = ap.offline_noise_var(:)';
            ap.backend = cfg.backend;
            ap.status = 'offline_profile_initialized';

            APProfile.band(b).ap(m) = ap;
        end
    end

    SpatialFP_init = SpatialFP;
    SpatialFP_init.cb_profile_initialized = true;
    SpatialFP_init.cb_profile_method = APProfile.method;
    SpatialFP_init = CB_recompute_fingerprint_derivatives_from_dBm(SpatialFP_init, Config);

    if cfg.verbose
        fprintf('[CB] APProfile initialized: B=%d, M=%d, G=%d, backend=%s\n', ...
            B, M, G, cfg.backend);
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
        P0 = NaN;
    end
end


function mean_model = fit_log_distance_mean(grid_xy, y0, ap_idx, SpatialFP, Config, band_id)
    ap_xy = [];
    if isfield(SpatialFP, 'APs') && isfield(SpatialFP.APs, 'pos_xy')
        ap_xy = SpatialFP.APs.pos_xy(ap_idx, :);
    elseif isfield(SpatialFP, 'ap_pos_xy')
        ap_xy = SpatialFP.ap_pos_xy(ap_idx, :);
    elseif isfield(Config, 'ap') && isfield(Config.ap, 'pos_xy') && ...
            size(Config.ap.pos_xy, 1) >= ap_idx
        ap_xy = Config.ap.pos_xy(ap_idx, :);
    elseif isfield(Config, 'calib_ap_profile') && ...
            isfield(Config.calib_ap_profile, 'ap_pos_xy') && ...
            size(Config.calib_ap_profile.ap_pos_xy, 1) >= ap_idx
        ap_xy = Config.calib_ap_profile.ap_pos_xy(ap_idx, :);
    end

    if isempty(ap_xy)
        ap_xy = [NaN, NaN];
        xfeat = zeros(size(y0));
        status = 'constant_mean_no_ap_xy';
    else
        d = sqrt(sum((grid_xy - ap_xy).^2, 2));
        d = max(d, 1.0);
        xfeat = log10(d);
        status = 'log_distance_mean';
    end

    valid = isfinite(y0) & isfinite(xfeat);
    if sum(valid) >= 2
        X = [ones(sum(valid), 1), xfeat(valid)];
        beta = X \ y0(valid);
    else
        beta = [nanmean_local(y0), 0];
        status = 'constant_mean_fallback';
    end

    mean_model = struct();
    mean_model.type = status;
    mean_model.ap_idx = ap_idx;
    mean_model.band_id = band_id;
    mean_model.ap_xy = ap_xy;
    mean_model.intercept_dB = beta(1);
    mean_model.slope_dB_per_log10m = beta(2);
    mean_model.n_fit = -beta(2) / 10;
end


function kp = make_kernel_params(cfg, b, m)
    kp = struct();
    kp.sigma_f_dB = get_param_value(cfg.sigma_f_dB, b, m);
    kp.ell_x_m = get_param_value(cfg.ell_x_m, b, m);
    kp.ell_y_m = get_param_value(cfg.ell_y_m, b, m);
    kp.eps_dist = 1e-6;
    kp.trajectory_quad_n = cfg.trajectory_quad_n;
    kp.rect_quad_n = cfg.rect_quad_n;
    kp.gaussian_quad_n = cfg.gaussian_quad_n;
end


function [kp_best, fit] = calibrate_offline_profile_hyperparams(grid_xy, y0, mean_model, kp0, cfg, b, m)
    kp_best = kp0;
    fit = struct();
    fit.enable = cfg.offline_fit_enable;
    fit.status = 'disabled';
    fit.band_id = b;
    fit.ap_id = m;
    fit.sigma0_dB = get_param_value(cfg.sigma0_dB, b, m);
    fit.best_fitness = NaN;
    fit.best_rmse_dB = NaN;
    fit.best_loglik_per_sample = NaN;
    fit.n_train = 0;
    fit.n_val = 0;

    valid = find(isfinite(y0) & all(isfinite(grid_xy), 2));
    if ~cfg.offline_fit_enable
        return;
    end
    if numel(valid) < cfg.offline_fit_min_points
        fit.status = 'not_enough_points';
        return;
    end

    old_rng = rng;
    rng(cfg.offline_fit_seed + 1000 * b + m);
    valid = valid(randperm(numel(valid)));
    rng(old_rng);

    n_train_full = max(2, round(cfg.offline_fit_train_fraction * numel(valid)));
    n_train = min(n_train_full, cfg.offline_fit_max_train_points);
    n_val = min(numel(valid) - n_train_full, cfg.offline_fit_max_val_points);
    if n_val < 2
        n_val = min(numel(valid) - n_train, cfg.offline_fit_max_val_points);
    end
    if n_train < 2 || n_val < 2
        fit.status = 'not_enough_split_points';
        return;
    end

    train_idx = valid(1:n_train);
    val_start = min(n_train_full + 1, numel(valid) - n_val + 1);
    val_idx = valid(val_start:(val_start + n_val - 1));

    Xtr = grid_xy(train_idx, :);
    Ytr = y0(train_idx);
    Xva = grid_xy(val_idx, :);
    Yva = y0(val_idx);
    mean_model_fit = fit_log_distance_mean_from_xy(Xtr, Ytr, mean_model, b, m);
    Mtr = CB_eval_ap_profile_mean(mean_model_fit, Xtr);
    Mva = CB_eval_ap_profile_mean(mean_model_fit, Xva);
    Rtr = Ytr(:) - Mtr(:);

    ell_grid = cfg.offline_fit_ell_grid_m(:)';
    sigma_f_grid = cfg.offline_fit_sigma_f_grid_dB(:)';
    sigma0_grid = cfg.offline_fit_sigma0_grid_dB(:)';

    best = struct('fitness', -inf, 'rmse', inf, 'loglik', -inf, ...
        'ell', kp0.ell_x_m, 'sigma_f', kp0.sigma_f_dB, 'sigma0', fit.sigma0_dB);

    for ie = 1:numel(ell_grid)
        ell = ell_grid(ie);
        for isf = 1:numel(sigma_f_grid)
            sigma_f = sigma_f_grid(isf);
            K0 = se_cov_xy(Xtr, Xtr, ell, ell, sigma_f^2);
            Kva = se_cov_xy(Xva, Xtr, ell, ell, sigma_f^2);
            for isn = 1:numel(sigma0_grid)
                sigma0 = sigma0_grid(isn);
                K = K0 + (sigma0^2 + cfg.jitter_abs) * eye(n_train);
                [L, ok] = stable_chol_offline(K);
                if ~ok
                    continue;
                end
                alpha = L' \ (L \ Rtr);
                pred = Mva + Kva * alpha;
                rmse = sqrt(mean((Yva(:) - pred(:)).^2));
                loglik = -0.5 * (Rtr' * alpha) - sum(log(diag(L))) - 0.5 * n_train * log(2*pi);
                loglik_per = loglik / n_train;
                fitness = cfg.offline_fit_w_rmse / max(rmse, eps) + ...
                    cfg.offline_fit_w_loglik * loglik_per;
                if fitness > best.fitness
                    best.fitness = fitness;
                    best.rmse = rmse;
                    best.loglik = loglik_per;
                    best.ell = ell;
                    best.sigma_f = sigma_f;
                    best.sigma0 = sigma0;
                end
            end
        end
    end

    if ~isfinite(best.fitness)
        fit.status = 'all_candidates_failed';
        return;
    end

    kp_best.sigma_f_dB = best.sigma_f;
    kp_best.ell_x_m = best.ell;
    kp_best.ell_y_m = best.ell;

    fit.status = 'ok';
    fit.train_fraction = cfg.offline_fit_train_fraction;
    fit.n_train = n_train;
    fit.n_val = n_val;
    fit.sigma_f_dB = best.sigma_f;
    fit.ell_x_m = best.ell;
    fit.ell_y_m = best.ell;
    fit.sigma0_dB = best.sigma0;
    fit.best_fitness = best.fitness;
    fit.best_rmse_dB = best.rmse;
    fit.best_loglik_per_sample = best.loglik;
    fit.mean_model_used_for_fitness = mean_model_fit;

    if cfg.offline_fit_verbose
        fprintf('[CB-offline-fit] band %d AP %d: ell=%.2f sigma_f=%.2f sigma0=%.2f rmse=%.3f fit=%.4f\n', ...
            b, m, best.ell, best.sigma_f, best.sigma0, best.rmse, best.fitness);
    end
end


function mean_model = fit_log_distance_mean_from_xy(xy, y, template, b, m)
    mean_model = template;
    mean_model.band_id = b;
    mean_model.ap_idx = m;
    if isfield(template, 'ap_xy') && all(isfinite(template.ap_xy))
        d = sqrt(sum((xy - template.ap_xy).^2, 2));
        d = max(d, 1.0);
        xfeat = log10(d);
        status = 'log_distance_mean_train_split';
    else
        xfeat = zeros(size(y));
        status = 'constant_mean_train_split';
    end

    valid = isfinite(y) & isfinite(xfeat);
    if sum(valid) >= 2
        X = [ones(sum(valid), 1), xfeat(valid)];
        beta = X \ y(valid);
    else
        beta = [nanmean_local(y), 0];
    end
    mean_model.type = status;
    mean_model.intercept_dB = beta(1);
    mean_model.slope_dB_per_log10m = beta(2);
    mean_model.n_fit = -beta(2) / 10;
end


function K = se_cov_xy(X1, X2, ell_x, ell_y, sigma_f2)
    dx = bsxfun(@minus, X1(:,1) / ell_x, (X2(:,1) / ell_x)');
    dy = bsxfun(@minus, X1(:,2) / ell_y, (X2(:,2) / ell_y)');
    K = sigma_f2 * exp(-0.5 * (dx.^2 + dy.^2));
end


function [L, ok] = stable_chol_offline(K)
    ok = false;
    K = 0.5 * (K + K');
    jitter = max(1e-9, 1e-8 * mean(max(diag(K), eps)));
    for i = 1:8
        [L, p] = chol(K + jitter * eye(size(K, 1)), 'lower');
        if p == 0
            ok = true;
            return;
        end
        jitter = jitter * 10;
    end
    L = [];
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


function v = nanmean_local(x)
    x = x(isfinite(x));
    if isempty(x)
        v = 0;
    else
        v = mean(x);
    end
end
