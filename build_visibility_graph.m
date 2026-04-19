function G = build_visibility_graph(src, dst, buildings)
% BUILD_VISIBILITY_GRAPH  节点 = {src, dst, 所有建筑角点(带 inset)}
%   边：两节点连线"不穿入任何矩形内部"则加边，权重=欧氏距离
%
%   为避免节点贴在矩形边界上与判定冲突，角点向外缩 inset

    inset = 0.05;  % m
    src = src(:)';
    dst = dst(:)';

    % 收集角点（外推）
    N = numel(buildings);
    corners = zeros(0, 2);
    for k = 1:N
        b = buildings(k);
        corners(end+1,:) = [b.xmin - inset, b.ymin - inset]; %#ok<AGROW>
        corners(end+1,:) = [b.xmax + inset, b.ymin - inset]; %#ok<AGROW>
        corners(end+1,:) = [b.xmax + inset, b.ymax + inset]; %#ok<AGROW>
        corners(end+1,:) = [b.xmin - inset, b.ymax + inset]; %#ok<AGROW>
    end

    nodes = [src; dst; corners];
    Nn = size(nodes, 1);

    % 邻接矩阵（上三角展开）
    A = inf(Nn, Nn);

    for i = 1:Nn
        A(i, i) = 0;
        for j = i+1:Nn
            [ok, len] = segment_passes_buildings(nodes(i,:), nodes(j,:), buildings);
            if ok
                A(i, j) = len;
                A(j, i) = len;
            end
        end
    end

    G.nodes = nodes;
    G.A = A;
    G.src_idx = 1;
    G.dst_idx = 2;
end


function [ok, len] = segment_passes_buildings(p1, p2, buildings)
% 若线段不"穿入"任何建筑内部，返回 ok=true
    len = norm(p2 - p1);
    if len < 1e-9
        ok = true;
        return;
    end
    ok = true;
    for k = 1:numel(buildings)
        if does_segment_enter_rect_local(p1, p2, buildings(k))
            % 允许：起点或终点本就在该建筑内 → 允许（否则终点在建筑内无路可走）
            if point_in_rect_local(p1, buildings(k)) || point_in_rect_local(p2, buildings(k))
                continue;
            end
            ok = false;
            return;
        end
    end
end


function tf = does_segment_enter_rect_local(p1, p2, bld)
    eps_geom = 1e-6;
    ts = linspace(0.02, 0.98, 13);
    tf = false;
    for i = 1:numel(ts)
        q = p1 + ts(i) * (p2 - p1);
        if (q(1) > bld.xmin + eps_geom) && (q(1) < bld.xmax - eps_geom) && ...
           (q(2) > bld.ymin + eps_geom) && (q(2) < bld.ymax - eps_geom)
            tf = true;
            return;
        end
    end
end


function tf = point_in_rect_local(p, bld)
    eps_geom = 1e-3;
    tf = (p(1) >= bld.xmin - eps_geom) && (p(1) <= bld.xmax + eps_geom) && ...
         (p(2) >= bld.ymin - eps_geom) && (p(2) <= bld.ymax + eps_geom);
end
