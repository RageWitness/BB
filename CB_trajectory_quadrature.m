function [pts, w] = CB_trajectory_quadrature(q, n_default)
% CB_TRAJECTORY_QUADRATURE  Deterministic quadrature points on a polyline.

    if nargin < 2 || isempty(n_default)
        n_default = 25;
    end

    if ~isfield(q, 'path_points') || size(q.path_points, 1) < 2
        error('[CB] trajectory prior requires path_points K x 2 with K>=2');
    end

    path = q.path_points;
    seg = diff(path, 1, 1);
    seg_len = sqrt(sum(seg.^2, 2));
    L = sum(seg_len);
    if L <= 0
        error('[CB] trajectory path length is zero');
    end

    n = n_default;
    if isfield(q, 'quad_n') && ~isempty(q.quad_n)
        n = q.quad_n;
    end
    n = max(3, n);

    s = linspace(0, L, n);
    pts = zeros(n, 2);
    cum_len = [0; cumsum(seg_len(:))];
    for i = 1:n
        si = s(i);
        k = find(cum_len <= si, 1, 'last');
        k = min(k, numel(seg_len));
        if seg_len(k) <= 0
            tau = 0;
        else
            tau = (si - cum_len(k)) / seg_len(k);
        end
        pts(i, :) = path(k, :) + tau * seg(k, :);
    end

    switch get_field_default(q, 'pdf_type', 'uniform_arclength')
        case 'uniform_arclength'
            w = ones(n, 1);
        case 'sampled_pdf'
            pdf_vals = get_field_default(q, 'pdf_values', []);
            if isempty(pdf_vals) || numel(pdf_vals) ~= n
                w = ones(n, 1);
            else
                w = pdf_vals(:);
            end
        case 'gaussian_on_s'
            params = get_field_default(q, 'pdf_params', struct());
            mu_s = get_field_default(params, 'mu_s', 0.5 * L);
            sig_s = get_field_default(params, 'sigma_s', 0.2 * L);
            w = exp(-0.5 * ((s(:) - mu_s) ./ max(sig_s, eps)).^2);
        otherwise
            w = ones(n, 1);
    end

    w(~isfinite(w) | w < 0) = 0;
    if sum(w) <= 0
        w = ones(n, 1);
    end
    w = w / sum(w);
end


function v = get_field_default(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end
