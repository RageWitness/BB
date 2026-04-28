function cfg = CD_uncertain_rss_defaults(Config)
% CD_UNCERTAIN_RSS_DEFAULTS  Defaults for opportunistic uncertain RSS calibration.

    if nargin < 1 || isempty(Config)
        Config = struct();
    end

    cfg = struct();
    if isfield(Config, 'cd') && isstruct(Config.cd)
        cfg = Config.cd;
    end

    cfg = set_default(cfg, 'enable', true);
    cfg = set_default(cfg, 'source_filter', 'opportunistic');
    cfg = set_default(cfg, 'supported_prior_types', {'region', 'gaussian'});

    cfg = set_default(cfg, 'patch_radius_m', 80);
    cfg = set_default(cfg, 'patch_taper_width_m', 20);
    cfg = set_default(cfg, 'max_base_points', 8);
    cfg = set_default(cfg, 'max_ray_points', 160);
    cfg = set_default(cfg, 'ray_width_m', 9);

    cfg = set_default(cfg, 'n_min', 0.5);
    cfg = set_default(cfg, 'n_max', 6.0);
    cfg = set_default(cfg, 'PL0_min_dB', -20);
    cfg = set_default(cfg, 'PL0_max_dB', 120);
    cfg = set_default(cfg, 'sigma_old_dB', 4.0);
    cfg = set_default(cfg, 'sigma_new_dB', 0.8);
    cfg = set_default(cfg, 'sigma_f_min_dB', 0.1);
    cfg = set_default(cfg, 'sigma_f_max_dB', 20);
    cfg = set_default(cfg, 'ell_min_m', 2);
    cfg = set_default(cfg, 'ell_max_m', 120);
    cfg = set_default(cfg, 'init_sigma_f_dB', 2.0);
    cfg = set_default(cfg, 'init_ell_m', 15);
    cfg = set_default(cfg, 'n0_mode', 'nearest_ap_per_band');

    cfg = set_default(cfg, 'blend_old_weight', 1.0);
    cfg = set_default(cfg, 'calibration_weight', 2.0);
    cfg = set_default(cfg, 'max_delta_dB', 25);
    cfg = set_default(cfg, 'min_ray_length_m', 1.0);
    cfg = set_default(cfg, 'gaussian_conf_chi2', 5.991);
    cfg = set_default(cfg, 'eps_dist', 1.0);
    cfg = set_default(cfg, 'verbose', true);
end


function s = set_default(s, name, value)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = value;
    end
end
