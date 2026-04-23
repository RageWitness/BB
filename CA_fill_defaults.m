function Config = CA_fill_defaults(Config)
% CA_FILL_DEFAULTS  补齐 Config.ca 所有默认参数（不覆盖已设字段）

    if ~isfield(Config, 'ca')
        Config.ca = struct();
    end
    ca = Config.ca;

    % --- interp ---
    ca = ensure_sub(ca, 'interp', struct());
    ca.interp = d(ca.interp, 'mode',  'idw');
    ca.interp = d(ca.interp, 'k',     4);
    ca.interp = d(ca.interp, 'power', 2);
    ca.interp = d(ca.interp, 'eps',   1e-9);

    % --- gpr ---
    ca = ensure_sub(ca, 'gpr', struct());
    ca.gpr = d(ca.gpr, 'n_min',             3);
    ca.gpr = d(ca.gpr, 'n_hyper_min',        6);
    ca.gpr = d(ca.gpr, 'jitter_init',        1e-6);
    ca.gpr = d(ca.gpr, 'jitter_max',         1e-2);
    ca.gpr = d(ca.gpr, 'jitter_factor',      10);
    ca.gpr = d(ca.gpr, 'max_opt_iter',       200);
    ca.gpr = d(ca.gpr, 'learn_hyperparams',  true);

    ca.gpr = ensure_sub(ca.gpr, 'init', struct());
    ca.gpr.init = d(ca.gpr.init, 'length_scale_x', 30);
    ca.gpr.init = d(ca.gpr.init, 'length_scale_y', 30);
    ca.gpr.init = d(ca.gpr.init, 'sigma_f',         2);
    ca.gpr.init = d(ca.gpr.init, 'sigma_n',         1);

    ca.gpr = ensure_sub(ca.gpr, 'bounds', struct());
    ca.gpr.bounds = d(ca.gpr.bounds, 'length_scale_min', 5);
    ca.gpr.bounds = d(ca.gpr.bounds, 'length_scale_max', 200);
    ca.gpr.bounds = d(ca.gpr.bounds, 'sigma_f_min',      0.1);
    ca.gpr.bounds = d(ca.gpr.bounds, 'sigma_f_max',      20);
    ca.gpr.bounds = d(ca.gpr.bounds, 'sigma_n_min',      0.1);
    ca.gpr.bounds = d(ca.gpr.bounds, 'sigma_n_max',      10);

    % --- kmeans ---
    ca = ensure_sub(ca, 'kmeans', struct());
    ca.kmeans = d(ca.kmeans, 'K',             4);
    ca.kmeans = d(ca.kmeans, 'max_iter',      100);
    ca.kmeans = d(ca.kmeans, 'n_replicates',  5);
    ca.kmeans = d(ca.kmeans, 'enable_blend',  true);
    ca.kmeans = d(ca.kmeans, 'blend_k',       2);
    ca.kmeans = d(ca.kmeans, 'blend_radius',  50);

    % --- update ---
    ca = ensure_sub(ca, 'update', struct());
    ca.update = d(ca.update, 'lambda_min',       0);
    ca.update = d(ca.update, 'lambda_max',       0.7);
    ca.update = d(ca.update, 'sigma_ref_dB',     3);
    ca.update = d(ca.update, 'n_ref',            5);
    ca.update = d(ca.update, 'max_delta_dB',     15);
    ca.update = d(ca.update, 'enable_adaptive_lambda', true);

    % --- sample ---
    ca = ensure_sub(ca, 'sample', struct());
    ca.sample = d(ca.sample, 'min_valid_aps',  2);
    ca.sample = d(ca.sample, 'skip_var',       1e6);

    % --- weight (frame-based CA) ---
    ca = ensure_sub(ca, 'weight', struct());
    ca.weight = d(ca.weight, 'tau_rec',    80);
    ca.weight = d(ca.weight, 'eps_w',      1e-6);

    % --- gpr n_eff thresholds (frame-based) ---
    ca.gpr = d(ca.gpr, 'n_min_eff',         3);
    ca.gpr = d(ca.gpr, 'n_hyper_min_eff',   6);

    % --- perturbation (地图老化扰动) ---
    ca = ensure_sub(ca, 'perturb', struct());
    ca.perturb = d(ca.perturb, 'enable',          false);
    ca.perturb = d(ca.perturb, 'mode',            'gp');
    ca.perturb = d(ca.perturb, 'amplitude_dB',    5);
    ca.perturb = d(ca.perturb, 'length_scale_m',  40);
    ca.perturb = d(ca.perturb, 'per_ap_std_dB',   2);
    ca.perturb = d(ca.perturb, 'seed',            123);

    Config.ca = ca;
end


function s = d(s, fname, val)
    if ~isfield(s, fname)
        s.(fname) = val;
    end
end

function s = ensure_sub(s, fname, val)
    if ~isfield(s, fname)
        s.(fname) = val;
    end
end
