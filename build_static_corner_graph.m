function VG = build_static_corner_graph(buildings)
% BUILD_STATIC_CORNER_GRAPH  预构建建筑角点的可见图（与 src/dst 无关）
%
%   输出：
%     VG.corners      - (Nc x 2) 角点坐标（向外缩 inset）
%     VG.A_cc         - (Nc x Nc) 角点-角点可见距离矩阵（不可见为 inf）
%     VG.buildings    - 建筑列表的拷贝
%     VG.inset        - 角点外推距离

    inset = 0.05;
    N = numel(buildings);
    corners = zeros(4*N, 2);
    for k = 1:N
        b = buildings(k);
        idx = 4*(k-1);
        corners(idx+1, :) = [b.xmin - inset, b.ymin - inset];
        corners(idx+2, :) = [b.xmax + inset, b.ymin - inset];
        corners(idx+3, :) = [b.xmax + inset, b.ymax + inset];
        corners(idx+4, :) = [b.xmin - inset, b.ymax + inset];
    end
    Nc = size(corners, 1);

    A_cc = inf(Nc, Nc);
    for i = 1:Nc
        A_cc(i, i) = 0;
        for j = i+1:Nc
            [ok, len] = segment_passes_buildings_static(corners(i,:), corners(j,:), buildings);
            if ok
                A_cc(i, j) = len;
                A_cc(j, i) = len;
            end
        end
    end

    VG.corners   = corners;
    VG.A_cc      = A_cc;
    VG.buildings = buildings;
    VG.inset     = inset;
    VG.Nc        = Nc;
end


function [ok, len] = segment_passes_buildings_static(p1, p2, buildings)
    len = norm(p2 - p1);
    ok = true;
    if len < 1e-9, return; end
    for k = 1:numel(buildings)
        if seg_enters_rect(p1, p2, buildings(k))
            if pt_on_or_in_rect(p1, buildings(k)) || pt_on_or_in_rect(p2, buildings(k))
                continue;
            end
            ok = false;
            return;
        end
    end
end

function tf = seg_enters_rect(p1, p2, bld)
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

function tf = pt_on_or_in_rect(p, bld)
    eps_geom = 1e-3;
    tf = (p(1) >= bld.xmin - eps_geom) && (p(1) <= bld.xmax + eps_geom) && ...
         (p(2) >= bld.ymin - eps_geom) && (p(2) <= bld.ymax + eps_geom);
end
