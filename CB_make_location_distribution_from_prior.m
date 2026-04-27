function q = CB_make_location_distribution_from_prior(location_prior)
% CB_MAKE_LOCATION_DISTRIBUTION_FROM_PRIOR  Convert source prior to q struct.

    if ~isstruct(location_prior) || ~isfield(location_prior, 'type')
        error('[CB] location_prior must be a struct with field type');
    end

    switch lower(location_prior.type)
        case 'exact'
            xy = get_prior_xy(location_prior);
            q = struct('type', 'point', 'x', xy(:)');

        case 'region'
            bbox = get_prior_bbox(location_prior);
            q = struct('type', 'rect', 'bbox', normalize_bbox(bbox));

        case 'gaussian'
            [mu, Sigma] = get_prior_gaussian(location_prior);
            q = struct('type', 'gaussian', 'mu', mu(:)', 'Sigma', Sigma);

        case 'trajectory'
            q = get_prior_trajectory(location_prior);

        otherwise
            error('[CB] unsupported location_prior.type=%s', location_prior.type);
    end
end


function xy = get_prior_xy(lp)
    if isfield(lp, 'xy') && numel(lp.xy) == 2
        xy = lp.xy;
    elseif isfield(lp, 'value') && numel(lp.value) == 2
        xy = lp.value;
    else
        error('[CB] exact location prior requires xy or value [x,y]');
    end
end


function bbox = get_prior_bbox(lp)
    v = [];
    if isfield(lp, 'bbox')
        v = lp.bbox;
    elseif isfield(lp, 'value')
        if isnumeric(lp.value)
            v = lp.value;
        elseif isstruct(lp.value)
            if isfield(lp.value, 'bbox')
                v = lp.value.bbox;
            elseif isfield(lp.value, 'region_geometry')
                v = lp.value.region_geometry;
            end
        end
    end
    if isempty(v) || numel(v) ~= 4
        error('[CB] region prior requires bbox [xmin,xmax,ymin,ymax]');
    end
    bbox = v(:)';
end


function bbox = normalize_bbox(bbox)
    x1 = min(bbox(1), bbox(2));
    x2 = max(bbox(1), bbox(2));
    y1 = min(bbox(3), bbox(4));
    y2 = max(bbox(3), bbox(4));
    bbox = [x1, x2, y1, y2];
end


function [mu, Sigma] = get_prior_gaussian(lp)
    if isfield(lp, 'mu')
        mu = lp.mu;
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'mu')
        mu = lp.value.mu;
    else
        error('[CB] gaussian prior requires mu');
    end

    if isfield(lp, 'Sigma')
        Sigma = lp.Sigma;
    elseif isfield(lp, 'sigma')
        Sigma = eye(2) * lp.sigma^2;
    elseif isfield(lp, 'value') && isstruct(lp.value)
        if isfield(lp.value, 'Sigma')
            Sigma = lp.value.Sigma;
        elseif isfield(lp.value, 'sigma')
            Sigma = eye(2) * lp.value.sigma^2;
        else
            error('[CB] gaussian prior requires Sigma or sigma');
        end
    else
        error('[CB] gaussian prior requires Sigma or sigma');
    end

    if isscalar(Sigma)
        Sigma = eye(2) * Sigma;
    elseif isvector(Sigma)
        Sigma = diag(Sigma(:));
    end
    Sigma = 0.5 * (Sigma + Sigma');
end


function q = get_prior_trajectory(lp)
    if isfield(lp, 'path_points')
        path = lp.path_points;
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'path_points')
        path = lp.value.path_points;
    else
        error('[CB] trajectory prior requires path_points');
    end
    q = struct();
    q.type = 'trajectory';
    q.path_points = path;
    q.pdf_type = get_nested_field(lp, 'pdf_type', 'uniform_arclength');
    q.pdf_params = get_nested_field(lp, 'pdf_params', struct());
    if isfield(lp, 'quad_n')
        q.quad_n = lp.quad_n;
    end
end


function v = get_nested_field(s, name, default_value)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    elseif isfield(s, 'value') && isstruct(s.value) && isfield(s.value, name) && ~isempty(s.value.(name))
        v = s.value.(name);
    else
        v = default_value;
    end
end
