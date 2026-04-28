function ray = CD_project_prior_to_ray(q, ap_xy, cfg)
% CD_PROJECT_PRIOR_TO_RAY  Project a 2D location prior onto one AP ray.

    ray = struct();
    ray.valid = false;
    ray.type = '';
    ray.reason = '';
    ray.ap_xy = ap_xy(:)';
    ray.dir = [NaN, NaN];
    ray.center_xy = [NaN, NaN];
    ray.mu_r = NaN;
    ray.sigma_r = NaN;
    ray.r_min = NaN;
    ray.r_max = NaN;

    if ~isfield(q, 'type')
        ray.reason = 'missing_prior_type';
        return;
    end

    switch lower(q.type)
        case 'gaussian'
            ray = project_gaussian(q, ray, cfg);
        case 'rect'
            ray = project_rect(q, ray, cfg);
        otherwise
            ray.reason = sprintf('unsupported_prior_%s', q.type);
    end
end


function ray = project_gaussian(q, ray, cfg)
    mu = q.mu(:)';
    z = mu - ray.ap_xy;
    r0 = norm(z);
    if r0 < cfg.min_ray_length_m
        ray.reason = 'gaussian_center_too_close_to_ap';
        return;
    end

    v = z / r0;
    S = 0.5 * (q.Sigma + q.Sigma');
    W = pinv(S);
    den = v * W * v';
    if den <= 0 || ~isfinite(den)
        ray.reason = 'invalid_gaussian_precision';
        return;
    end

    var_r = 1 / den;
    mu_r = var_r * (v * W * z');
    sigma_r = sqrt(max(var_r, 0));

    ray.valid = isfinite(mu_r) && mu_r > cfg.eps_dist;
    ray.type = 'gaussian';
    ray.reason = '';
    ray.dir = v;
    ray.center_xy = mu;
    ray.mu_r = max(mu_r, cfg.eps_dist);
    ray.sigma_r = sigma_r;
    ray.r_min = max(cfg.eps_dist, ray.mu_r - sqrt(cfg.gaussian_conf_chi2) * sigma_r);
    ray.r_max = max(ray.r_min, ray.mu_r + sqrt(cfg.gaussian_conf_chi2) * sigma_r);
    if ~ray.valid
        ray.reason = 'invalid_gaussian_projection';
    end
end


function ray = project_rect(q, ray, cfg)
    bbox = q.bbox;
    c = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
    z = c - ray.ap_xy;
    r0 = norm(z);
    if r0 < cfg.min_ray_length_m
        ray.reason = 'rect_center_too_close_to_ap';
        return;
    end
    v = z / r0;

    [r_min, r_max, ok] = intersect_ray_rect(ray.ap_xy, v, bbox);
    if ~ok || r_max <= cfg.eps_dist
        ray.reason = 'ray_misses_rect';
        return;
    end

    r_min = max(r_min, cfg.eps_dist);
    r_max = max(r_max, r_min);
    ray.valid = true;
    ray.type = 'rect';
    ray.reason = '';
    ray.dir = v;
    ray.center_xy = c;
    ray.r_min = r_min;
    ray.r_max = r_max;
    ray.mu_r = 0.5 * (r_min + r_max);
    ray.sigma_r = (r_max - r_min) / sqrt(12);
end


function [r_min, r_max, ok] = intersect_ray_rect(origin, dirv, bbox)
    xmin = bbox(1); xmax = bbox(2);
    ymin = bbox(3); ymax = bbox(4);

    t0 = -inf;
    t1 = inf;

    if abs(dirv(1)) < eps
        if origin(1) < xmin || origin(1) > xmax
            ok = false; r_min = NaN; r_max = NaN; return;
        end
    else
        tx1 = (xmin - origin(1)) / dirv(1);
        tx2 = (xmax - origin(1)) / dirv(1);
        t0 = max(t0, min(tx1, tx2));
        t1 = min(t1, max(tx1, tx2));
    end

    if abs(dirv(2)) < eps
        if origin(2) < ymin || origin(2) > ymax
            ok = false; r_min = NaN; r_max = NaN; return;
        end
    else
        ty1 = (ymin - origin(2)) / dirv(2);
        ty2 = (ymax - origin(2)) / dirv(2);
        t0 = max(t0, min(ty1, ty2));
        t1 = min(t1, max(ty1, ty2));
    end

    r_min = max(t0, 0);
    r_max = t1;
    ok = isfinite(r_min) && isfinite(r_max) && r_max >= r_min;
end
