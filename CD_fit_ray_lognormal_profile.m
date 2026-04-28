function fit = CD_fit_ray_lognormal_profile(r_old, y_old, P0_dBm, y_event_dBm, ray, cfg)
% CD_FIT_RAY_LOGNORMAL_PROFILE  Fit 1D RSS log-distance curve on one AP ray.
%
% Model in dBm:
%   y(r) = A + B*log10(r) + residual(r)
%   n0 = -B/10, PL0 = P0 - A

    fit = init_fit();

    valid = isfinite(r_old) & isfinite(y_old) & r_old > cfg.eps_dist;
    r_old = r_old(valid);
    y_old = y_old(valid);
    if isempty(r_old)
        fit.warning = 'no_old_ray_samples';
        return;
    end

    if numel(r_old) > cfg.max_ray_points
        idx = round(linspace(1, numel(r_old), cfg.max_ray_points));
        r_old = r_old(idx);
        y_old = y_old(idx);
    end

    A0 = P0_dBm - 20;
    B0 = -25;
    z0 = zeros(1, 5);
    z0(1) = to_unbounded(A0, P0_dBm - cfg.PL0_max_dB, P0_dBm - cfg.PL0_min_dB);
    z0(2) = to_unbounded(B0, -10 * cfg.n_max, -10 * cfg.n_min);
    if strcmp(ray.type, 'rect')
        z0(3) = to_unbounded(ray.mu_r, ray.r_min, ray.r_max);
    else
        z0(3) = to_unbounded(ray.mu_r, ray.r_min, ray.r_max);
    end
    z0(4) = to_unbounded(cfg.init_sigma_f_dB, cfg.sigma_f_min_dB, cfg.sigma_f_max_dB);
    z0(5) = to_unbounded(cfg.init_ell_m, cfg.ell_min_m, cfg.ell_max_m);

    opts = optimset('Display', 'off', 'MaxIter', 300, 'MaxFunEvals', 1200);
    obj = @(z) objective(z, r_old, y_old, P0_dBm, y_event_dBm, ray, cfg);
    [z_best, nmll] = fminsearch(obj, z0, opts);
    theta = decode_theta(z_best, P0_dBm, ray, cfg);

    if ~isfinite(nmll)
        fit.warning = 'optimization_failed';
        return;
    end

    fit.valid = true;
    fit.warning = '';
    fit.A_dBm = theta.A;
    fit.B_dB_per_log10m = theta.B;
    fit.n_hat = max(cfg.n_min, min(cfg.n_max, -theta.B / 10));
    fit.PL0_hat_dB = P0_dBm - theta.A;
    fit.r_event_hat = theta.r_event;
    fit.sigma_f_dB = theta.sigma_f;
    fit.ell_m = theta.ell;
    fit.nmll = nmll;
end


function fit = init_fit()
    fit = struct();
    fit.valid = false;
    fit.warning = 'not_run';
    fit.A_dBm = NaN;
    fit.B_dB_per_log10m = NaN;
    fit.n_hat = NaN;
    fit.PL0_hat_dB = NaN;
    fit.r_event_hat = NaN;
    fit.sigma_f_dB = NaN;
    fit.ell_m = NaN;
    fit.nmll = NaN;
end


function nmll = objective(z, r_old, y_old, P0_dBm, y_event_dBm, ray, cfg)
    th = decode_theta(z, P0_dBm, ray, cfg);

    r_all = [r_old(:); 1; th.r_event];
    y_all = [y_old(:); P0_dBm; y_event_dBm];
    mu = th.A + th.B * log10(max(r_all, cfg.eps_dist));
    res = y_all(:) - mu(:);

    D = abs(r_all(:) - r_all(:)');
    K = th.sigma_f^2 * exp(-0.5 * (D / max(th.ell, eps)).^2);

    n_old = numel(r_old);
    noise = cfg.sigma_old_dB^2 * ones(n_old + 2, 1);
    noise(n_old + 1) = cfg.sigma_new_dB^2;
    if strcmp(ray.type, 'gaussian')
        noise(n_old + 2) = cfg.sigma_new_dB^2 + ...
            (th.B / (max(th.r_event, cfg.eps_dist) * log(10)))^2 * ray.sigma_r^2;
    else
        noise(n_old + 2) = cfg.sigma_new_dB^2;
    end

    C = K + diag(noise);
    [L, ok] = stable_chol(C);
    if ~ok
        nmll = inf;
        return;
    end
    alpha = L' \ (L \ res);
    nmll = 0.5 * (res' * alpha) + sum(log(diag(L))) + 0.5 * numel(res) * log(2*pi);
end


function th = decode_theta(z, P0_dBm, ray, cfg)
    th = struct();
    th.A = from_unbounded(z(1), P0_dBm - cfg.PL0_max_dB, P0_dBm - cfg.PL0_min_dB);
    th.B = from_unbounded(z(2), -10 * cfg.n_max, -10 * cfg.n_min);
    if strcmp(ray.type, 'gaussian')
        th.r_event = max(ray.mu_r, cfg.eps_dist);
    else
        th.r_event = from_unbounded(z(3), ray.r_min, ray.r_max);
    end
    th.sigma_f = from_unbounded(z(4), cfg.sigma_f_min_dB, cfg.sigma_f_max_dB);
    th.ell = from_unbounded(z(5), cfg.ell_min_m, cfg.ell_max_m);
end


function z = to_unbounded(x, lo, hi)
    if hi <= lo + eps
        z = 0;
        return;
    end
    x = min(max(x, lo + eps), hi - eps);
    p = (x - lo) / (hi - lo);
    z = log(p / (1 - p));
end


function x = from_unbounded(z, lo, hi)
    if hi <= lo + eps
        x = 0.5 * (lo + hi);
        return;
    end
    p = 1 ./ (1 + exp(-z));
    x = lo + (hi - lo) * p;
end


function [L, ok] = stable_chol(C)
    ok = false;
    C = 0.5 * (C + C');
    jitter = max(1e-9, 1e-8 * mean(max(diag(C), eps)));
    for k = 1:8
        [L, p] = chol(C + jitter * eye(size(C, 1)), 'lower');
        if p == 0
            ok = true;
            return;
        end
        jitter = jitter * 10;
    end
    L = [];
end
