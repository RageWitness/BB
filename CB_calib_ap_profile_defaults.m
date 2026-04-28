function cfg = CB_calib_ap_profile_defaults(Config)
% CB_CALIB_AP_PROFILE_DEFAULTS  Defaults for AP-profile RSS calibration.

    if nargin < 1 || isempty(Config)
        Config = struct();
    end

    cfg = struct();
    if isfield(Config, 'calib_ap_profile') && isstruct(Config.calib_ap_profile)
        cfg = Config.calib_ap_profile;
    end

    cfg = set_default(cfg, 'enable', true);
    cfg = set_default(cfg, 'power_known_only', true);
    cfg = set_default(cfg, 'kernel', 'se_ard');
    cfg = set_default(cfg, 'ap_pos_xy', []);

    cfg = set_default(cfg, 'sigma_f_dB', 8);
    cfg = set_default(cfg, 'ell_x_m', 35);
    cfg = set_default(cfg, 'ell_y_m', 35);
    cfg = set_default(cfg, 'sigma0_dB', 8);
    cfg = set_default(cfg, 'sigma_meas_dB', 1.5);
    cfg = set_default(cfg, 'gamma_region_var', 1.0);
    cfg = set_default(cfg, 'gamma_gaussian_var', 0.0);

    cfg = set_default(cfg, 'trajectory_quad_n', 25);
    cfg = set_default(cfg, 'rect_quad_n', 9);
    cfg = set_default(cfg, 'gaussian_quad_n', 9);

    cfg = set_default(cfg, 'backend', 'scalable_subset_gp');
    cfg = set_default(cfg, 'max_train_points_exact_gp', 350);
    cfg = set_default(cfg, 'max_offline_train_points', 300);
    cfg = set_default(cfg, 'prediction_block_size', 1500);
    cfg = set_default(cfg, 'jitter_rel', 1e-6);
    cfg = set_default(cfg, 'jitter_abs', 1e-8);

    cfg = set_default(cfg, 'offline_fit_enable', true);
    cfg = set_default(cfg, 'offline_fit_train_fraction', 0.60);
    cfg = set_default(cfg, 'offline_fit_max_train_points', 180);
    cfg = set_default(cfg, 'offline_fit_max_val_points', 120);
    cfg = set_default(cfg, 'offline_fit_min_points', 30);
    cfg = set_default(cfg, 'offline_fit_seed', 20260428);
    cfg = set_default(cfg, 'offline_fit_ell_grid_m', [15, 30, 60, 100]);
    cfg = set_default(cfg, 'offline_fit_sigma_f_grid_dB', [2, 5, 8, 12]);
    cfg = set_default(cfg, 'offline_fit_sigma0_grid_dB', [0.8, 1.5, 3, 6, 10]);
    cfg = set_default(cfg, 'offline_fit_w_rmse', 1.0);
    cfg = set_default(cfg, 'offline_fit_w_loglik', 0.05);
    cfg = set_default(cfg, 'offline_fit_verbose', false);

    cfg = set_default(cfg, 'validation_fraction', 0.20);
    cfg = set_default(cfg, 'loc_rmse_threshold_m', 13);
    cfg = set_default(cfg, 'enable_fast_search', true);
    cfg = set_default(cfg, 'fast_search_backend', 'kernel_smoothing_fast');
    cfg = set_default(cfg, 'fast_ell_grid_m', [20, 40, 80, 120]);
    cfg = set_default(cfg, 'fast_trust_ratio_grid', [2, 4, 8, 12]);
    cfg = set_default(cfg, 'fast_refine_enable', true);
    cfg = set_default(cfg, 'fast_refine_scale', [0.75, 1.0, 1.25]);
    cfg = set_default(cfg, 'enable_ga_feedback', false);
    cfg = set_default(cfg, 'ga_population', 10);
    cfg = set_default(cfg, 'ga_generations', 5);
    cfg = set_default(cfg, 'ga_mutation_rate', 0.25);
    cfg = set_default(cfg, 'ga_crossover_rate', 0.65);
    cfg = set_default(cfg, 'ga_elite_count', 2);

    cfg = set_default(cfg, 'reject_power_none', true);
    cfg = set_default(cfg, 'verbose', true);
end


function s = set_default(s, name, value)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = value;
    end
end
