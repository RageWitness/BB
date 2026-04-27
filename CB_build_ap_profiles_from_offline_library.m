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

            ap = struct();
            ap.kernel_params = kp;
            ap.mean_model = mean_model;
            ap.offline_xy = grid_xy;
            ap.offline_y_dBm = y0;
            ap.offline_noise_var = get_param_value(cfg.sigma0_dB, b, m)^2 * ones(G, 1);
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
