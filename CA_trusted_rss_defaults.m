function cfg = CA_trusted_rss_defaults(Config)
% CA_TRUSTED_RSS_DEFAULTS  Defaults for trusted-source RF_raw calibration.

    cfg = struct();

    cfg.N_min = 1;
    cfg.N_max = 5;
    cfg.K_coarse = 100;
    cfg.TolX = 1e-6;
    cfg.eps_dist = 1e-3;
    cfg.eps_power = 1e-300;
    cfg.ill_cond_log_eta_thresh = 1e-3;
    cfg.power_equal_tol_dB = 1e-6;
    cfg.min_frames = 1;
    cfg.fusion_mode = 'mean';
    cfg.verbose = true;

    if nargin < 1 || isempty(Config)
        return;
    end

    if isfield(Config, 'ca') && isfield(Config.ca, 'trusted_rss')
        cfg = merge_struct(cfg, Config.ca.trusted_rss);
    end
end


function a = merge_struct(a, b)
    names = fieldnames(b);
    for i = 1:numel(names)
        a.(names{i}) = b.(names{i});
    end
end
