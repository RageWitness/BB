function result = CA_estimate_n0_top2_ratio_search(ap_xy, source_xy, rss_dB_frames, cfg)
% CA_ESTIMATE_N0_TOP2_RATIO_SEARCH  Estimate path-loss exponent from top-2 APs.
%
% Uses only the strongest and second strongest AP average RSS values.

    if nargin < 4 || isempty(cfg)
        cfg = CA_trusted_rss_defaults();
    else
        cfg = fill_missing_cfg(cfg);
    end

    result = init_result();

    ap_xy = double(ap_xy);
    source_xy = double(source_xy(:)');
    M = size(ap_xy, 1);

    if size(ap_xy, 2) ~= 2 || numel(source_xy) ~= 2
        result.warning = 'invalid_geometry_dimension';
        return;
    end

    [rss_mat, ok_shape] = normalize_rss_shape(rss_dB_frames, M);
    if ~ok_shape
        result.warning = 'invalid_rss_dimension';
        return;
    end

    Rbar_dB = nanmean_cols(rss_mat);
    Q_lin = 10.^(Rbar_dB / 10);
    result.Rbar_dB = Rbar_dB;
    result.Q_lin = Q_lin;

    valid_q = isfinite(Q_lin) & (Q_lin > cfg.eps_power);
    if sum(valid_q) < 2
        result.warning = 'fewer_than_two_valid_ap_powers';
        return;
    end

    Q_sort = Q_lin;
    Q_sort(~valid_q) = -inf;
    [~, idx_sorted] = sort(Q_sort, 'descend');
    r = idx_sorted(1);
    i = idx_sorted(2);

    d_all = sqrt(sum((ap_xy - source_xy).^2, 2));
    d_r = d_all(r);
    d_i = d_all(i);
    Q_r = Q_lin(r);
    Q_i = Q_lin(i);

    result.ref_ap = r;
    result.second_ap = i;

    if d_r <= cfg.eps_dist || d_i <= cfg.eps_dist
        result.warning = 'invalid_top2_distance';
        return;
    end
    if Q_r <= cfg.eps_power || Q_i <= cfg.eps_power
        result.warning = 'invalid_top2_power';
        return;
    end

    rho = d_r / d_i;
    eta = Q_i / Q_r;
    result.rho = rho;
    result.eta = eta;

    warn_list = {};
    if abs(log(eta)) < cfg.ill_cond_log_eta_thresh
        warn_list{end+1} = 'ill_conditioned_eta'; %#ok<AGROW>
    end

    u_min = 1 / cfg.N_max;
    u_max = 1 / cfg.N_min;
    if ~(isfinite(u_min) && isfinite(u_max) && u_min > 0 && u_max > u_min)
        result.warning = 'invalid_search_bounds';
        return;
    end

    K = max(3, round(cfg.K_coarse));
    u_grid = linspace(u_min, u_max, K);
    J_grid = (rho - eta .^ u_grid).^2;
    [~, j_best] = min(J_grid);

    j_lo = max(1, j_best - 1);
    j_hi = min(K, j_best + 1);
    lo = u_grid(j_lo);
    hi = u_grid(j_hi);
    if lo == hi
        lo = u_min;
        hi = u_max;
    end

    obj = @(u) (rho - eta.^u).^2;
    opts = optimset('TolX', cfg.TolX, 'Display', 'off');
    [u_hat, J_min] = fminbnd(obj, lo, hi, opts);
    n_hat = 1 / u_hat;

    if n_hat < cfg.N_min
        n_hat = cfg.N_min;
        u_hat = 1 / n_hat;
        warn_list{end+1} = 'n_hat_clipped_to_min'; %#ok<AGROW>
    elseif n_hat > cfg.N_max
        n_hat = cfg.N_max;
        u_hat = 1 / n_hat;
        warn_list{end+1} = 'n_hat_clipped_to_max'; %#ok<AGROW>
    end

    result.n_hat = n_hat;
    result.u_hat = u_hat;
    result.J_min = J_min;
    result.valid = true;
    result.warning = join_warnings(warn_list);
end


function cfg = fill_missing_cfg(cfg)
    d = CA_trusted_rss_defaults();
    names = fieldnames(d);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k}) || isempty(cfg.(names{k}))
            cfg.(names{k}) = d.(names{k});
        end
    end
end


function result = init_result()
    result = struct();
    result.n_hat = NaN;
    result.u_hat = NaN;
    result.ref_ap = NaN;
    result.second_ap = NaN;
    result.rho = NaN;
    result.eta = NaN;
    result.J_min = NaN;
    result.valid = false;
    result.warning = '';
    result.Rbar_dB = [];
    result.Q_lin = [];
end


function [rss_mat, ok] = normalize_rss_shape(rss_dB_frames, M)
    ok = true;
    rss = double(rss_dB_frames);
    if isvector(rss)
        if numel(rss) ~= M
            ok = false;
            rss_mat = [];
            return;
        end
        rss_mat = reshape(rss, 1, M);
        return;
    end

    if size(rss, 2) == M
        rss_mat = rss;
    elseif size(rss, 1) == M
        rss_mat = rss';
    else
        ok = false;
        rss_mat = [];
    end
end


function m = nanmean_cols(X)
    X(~isfinite(X)) = NaN;
    n = sum(~isnan(X), 1);
    s = sum(replace_nan(X, 0), 1);
    m = s ./ max(n, 1);
    m(n == 0) = NaN;
end


function X = replace_nan(X, val)
    X(isnan(X)) = val;
end


function s = join_warnings(warn_list)
    if isempty(warn_list)
        s = '';
    else
        s = warn_list{1};
        for k = 2:numel(warn_list)
            s = [s, ';', warn_list{k}]; %#ok<AGROW>
        end
    end
end
