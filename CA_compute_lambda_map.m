function DeltaA = CA_compute_lambda_map(DeltaA, Config)
% CA_COMPUTE_LAMBDA_MAP  根据 var 和 n_eff 计算自适应 lambda map
%
%   lambda(b,m,g) = lambda_max * c_var(b,m,g) * c_den(b,m,g)
%   clipped to [lambda_min, lambda_max]

    Config = CA_fill_defaults(Config);
    up = Config.ca.update;
    sigma_ref = up.sigma_ref_dB;
    n_ref     = up.n_ref;
    lam_min   = up.lambda_min;
    lam_max   = up.lambda_max;
    adaptive  = up.enable_adaptive_lambda;

    B = numel(DeltaA.band);
    for b = 1:B
        [M_ap, G] = size(DeltaA.band(b).mu);

        if ~adaptive
            lam = lam_max * ones(M_ap, G);
            lam(~DeltaA.band(b).valid) = 0;
            DeltaA.band(b).lambda = lam;
            continue;
        end

        var_map = DeltaA.band(b).var;
        neff    = DeltaA.band(b).n_eff;
        valid   = DeltaA.band(b).valid;

        c_var = sigma_ref^2 ./ (var_map + sigma_ref^2);
        c_den = min(1, neff / n_ref);

        lam = lam_max * c_var .* c_den;
        lam = max(lam_min, min(lam_max, lam));
        lam(~valid) = 0;

        DeltaA.band(b).lambda = lam;
    end
end
