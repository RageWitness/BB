function [theta, opt_info] = CA_optimize_gpr_hyperparams(X, y, Config)
% CA_OPTIMIZE_GPR_HYPERPARAMS  负对数边际似然优化 SE-ARD 核超参数
%
%   theta = [log(ell_x), log(ell_y), log(sigma_f), log(sigma_n)]
%   X: N x 2,  y: N x 1 (无 NaN)

    Config = CA_fill_defaults(Config);
    gpr_cfg = Config.ca.gpr;
    bnd = gpr_cfg.bounds;

    theta0 = [ log(gpr_cfg.init.length_scale_x), ...
               log(gpr_cfg.init.length_scale_y), ...
               log(gpr_cfg.init.sigma_f), ...
               log(gpr_cfg.init.sigma_n) ];

    lb = [ log(bnd.length_scale_min), log(bnd.length_scale_min), ...
           log(bnd.sigma_f_min),      log(bnd.sigma_n_min) ];
    ub = [ log(bnd.length_scale_max), log(bnd.length_scale_max), ...
           log(bnd.sigma_f_max),      log(bnd.sigma_n_max) ];

    jitter_init  = gpr_cfg.jitter_init;

    obj = @(th) nll_se_ard(th, X, y, jitter_init, gpr_cfg.jitter_max, gpr_cfg.jitter_factor, lb, ub);

    opts = optimset('MaxIter', gpr_cfg.max_opt_iter, 'MaxFunEvals', gpr_cfg.max_opt_iter * 4, ...
                    'Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-5);

    opt_info = struct('status', 'ok', 'nll_init', NaN, 'nll_final', NaN, 'n_iter', 0);

    try
        opt_info.nll_init = obj(theta0);
        [theta_raw, fval, exitflag, output] = fminsearch(obj, theta0, opts);
        opt_info.nll_final = fval;
        opt_info.n_iter = output.iterations;

        theta = max(lb, min(ub, theta_raw));

        if exitflag <= 0
            opt_info.status = 'fminsearch_did_not_converge';
        end
    catch ME
        opt_info.status = sprintf('optimization_failed: %s', ME.message);
        theta = theta0;
    end
end


function nll = nll_se_ard(theta, X, y, jitter, jitter_max, jitter_factor, lb, ub)
    theta_c = max(lb, min(ub, theta));

    ell_x = exp(theta_c(1));
    ell_y = exp(theta_c(2));
    sf    = exp(theta_c(3));
    sn    = exp(theta_c(4));

    N = size(X, 1);
    D = (X(:,1) - X(:,1)').^2 / ell_x^2 + (X(:,2) - X(:,2)').^2 / ell_y^2;
    K = sf^2 * exp(-0.5 * D) + sn^2 * eye(N);

    % soft penalty for bound violation
    penalty = 100 * sum((theta - theta_c).^2);

    j = jitter;
    while j <= jitter_max
        try
            L = chol(K + j * eye(N), 'lower');
            alpha = L' \ (L \ y);
            nll = 0.5 * (y' * alpha) + sum(log(diag(L))) + 0.5 * N * log(2*pi) + penalty;
            return;
        catch
            j = j * jitter_factor;
        end
    end

    nll = 1e10 + penalty;
end
