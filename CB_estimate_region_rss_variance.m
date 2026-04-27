function [var_dB2, info] = CB_estimate_region_rss_variance(grid_xy, R_init, q_rect, cfg)
% CB_ESTIMATE_REGION_RSS_VARIANCE  Estimate RSS variance inside a rectangle.

    bbox = q_rect.bbox(:)';
    inside = grid_xy(:,1) >= bbox(1) & grid_xy(:,1) <= bbox(2) & ...
             grid_xy(:,2) >= bbox(3) & grid_xy(:,2) <= bbox(4);

    vals = R_init(inside);
    if isempty(vals)
        center = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
        d2 = sum((grid_xy - center).^2, 2);
        k = min(numel(d2), max(10, cfg.rect_quad_n^2));
        [~, ord] = sort(d2, 'ascend');
        vals = R_init(ord(1:k));
    end

    vals = vals(isfinite(vals));
    if numel(vals) <= 1
        var_dB2 = 0;
    else
        var_dB2 = var(vals, 0);
    end

    info = struct();
    info.bbox = bbox;
    info.area = max(bbox(2)-bbox(1), 0) * max(bbox(4)-bbox(3), 0);
    info.n_grid_used = numel(vals);
    info.mean_R_dBm = mean(vals);
    info.var_R_dB2 = var_dB2;
end
