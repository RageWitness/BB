function pathInfo = find_diffraction_path(src, dst, buildings, VG)
% FIND_DIFFRACTION_PATH  Dijkstra 在 (src, dst, 角点) 图上找最短路径
%
%   传入 VG（由 build_static_corner_graph 预构建）以避免每次重算角点-角点边

    if nargin < 4 || isempty(VG)
        VG = build_static_corner_graph(buildings);
    end

    src = src(:)';
    dst = dst(:)';
    corners = VG.corners;
    Nc = VG.Nc;
    Nn = Nc + 2;

    % 节点：1=src, 2=dst, 3..Nn=corners
    A = inf(Nn, Nn);
    A(3:Nn, 3:Nn) = VG.A_cc;
    for i = 1:Nn, A(i, i) = 0; end

    % src/dst 对所有角点的可见边
    for i = 1:Nc
        [ok, len] = seg_pass(src, corners(i, :), buildings);
        if ok, A(1, 2+i) = len; A(2+i, 1) = len; end

        [ok, len] = seg_pass(dst, corners(i, :), buildings);
        if ok, A(2, 2+i) = len; A(2+i, 2) = len; end
    end

    % src - dst 直接边
    [ok, len] = seg_pass(src, dst, buildings);
    if ok, A(1,2) = len; A(2,1) = len; end

    % Dijkstra (src=1, dst=2)
    dist = inf(1, Nn);
    prev = zeros(1, Nn);
    visited = false(1, Nn);
    dist(1) = 0;

    for step = 1:Nn
        tmp = dist; tmp(visited) = inf;
        [du, u] = min(tmp);
        if isinf(du), break; end
        visited(u) = true;
        if u == 2, break; end
        cand = find(~visited & A(u, :) < inf);
        for v = cand
            alt = du + A(u, v);
            if alt < dist(v)
                dist(v) = alt; prev(v) = u;
            end
        end
    end

    pathInfo = struct('ok', false, 'path_points', [], ...
        'total_length', inf, 'num_corners', 0, ...
        'wall_crossings_along_path', 0);

    if isinf(dist(2)), return; end

    seq = 2; cur = 2;
    while prev(cur) ~= 0
        cur = prev(cur);
        seq(end+1) = cur; %#ok<AGROW>
    end
    seq = fliplr(seq);

    nodes = [src; dst; corners];
    pts = nodes(seq, :);
    pathInfo.path_points = pts;
    pathInfo.total_length = dist(2);
    pathInfo.num_corners = max(0, size(pts, 1) - 2);
    pathInfo.ok = true;

    w_total = 0;
    for i = 1:size(pts, 1) - 1
        w = count_wall_crossings(pts(i,:), pts(i+1,:), buildings);
        w_total = w_total + w;
    end
    pathInfo.wall_crossings_along_path = w_total;
end


function [ok, len] = seg_pass(p1, p2, buildings)
    len = norm(p2 - p1);
    ok = true;
    if len < 1e-9, return; end
    for k = 1:numel(buildings)
        if seg_enters(p1, p2, buildings(k))
            if pt_in_or_on(p1, buildings(k)) || pt_in_or_on(p2, buildings(k))
                continue;
            end
            ok = false; return;
        end
    end
end

function tf = seg_enters(p1, p2, bld)
    eps_geom = 1e-6;
    ts = linspace(0.02, 0.98, 13);
    tf = false;
    for i = 1:numel(ts)
        q = p1 + ts(i) * (p2 - p1);
        if (q(1) > bld.xmin + eps_geom) && (q(1) < bld.xmax - eps_geom) && ...
           (q(2) > bld.ymin + eps_geom) && (q(2) < bld.ymax - eps_geom)
            tf = true; return;
        end
    end
end

function tf = pt_in_or_on(p, bld)
    eps_geom = 1e-3;
    tf = (p(1) >= bld.xmin - eps_geom) && (p(1) <= bld.xmax + eps_geom) && ...
         (p(2) >= bld.ymin - eps_geom) && (p(2) <= bld.ymax + eps_geom);
end
