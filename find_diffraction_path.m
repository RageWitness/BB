function pathInfo = find_diffraction_path(src, dst, buildings)
% FIND_DIFFRACTION_PATH  角点可见图 + Dijkstra 最短路径
%
%   输出：
%     pathInfo.path_points              - (Np x 2) 路径点列（含起点终点）
%     pathInfo.total_length             - dc_i
%     pathInfo.num_corners              - Nc_i（中间节点的转向数）
%     pathInfo.wall_crossings_along_path- Nw_dif_i（沿途墙损累加）
%     pathInfo.ok                       - 是否找到路径

    G = build_visibility_graph(src, dst, buildings);
    A = G.A;
    Nn = size(A, 1);

    % Dijkstra
    dist = inf(1, Nn);
    prev = zeros(1, Nn);
    visited = false(1, Nn);
    dist(G.src_idx) = 0;

    for step = 1:Nn
        % 未访问中最小
        tmp = dist;
        tmp(visited) = inf;
        [du, u] = min(tmp);
        if isinf(du), break; end
        visited(u) = true;
        if u == G.dst_idx, break; end
        for v = 1:Nn
            if ~visited(v) && A(u, v) < inf
                alt = du + A(u, v);
                if alt < dist(v)
                    dist(v) = alt;
                    prev(v) = u;
                end
            end
        end
    end

    pathInfo = struct();
    pathInfo.ok = false;
    pathInfo.path_points = [];
    pathInfo.total_length = inf;
    pathInfo.num_corners = 0;
    pathInfo.wall_crossings_along_path = 0;

    if isinf(dist(G.dst_idx))
        return;
    end

    % 回溯
    seq = G.dst_idx;
    cur = G.dst_idx;
    while prev(cur) ~= 0
        cur = prev(cur);
        seq(end+1) = cur; %#ok<AGROW>
    end
    seq = fliplr(seq);

    pts = G.nodes(seq, :);
    pathInfo.path_points = pts;
    pathInfo.total_length = dist(G.dst_idx);
    pathInfo.num_corners = max(0, size(pts, 1) - 2);  % 去起点终点
    pathInfo.ok = true;

    % 沿途墙穿越累加
    w_total = 0;
    for i = 1:size(pts, 1) - 1
        w = count_wall_crossings(pts(i,:), pts(i+1,:), buildings);
        w_total = w_total + w;
    end
    pathInfo.wall_crossings_along_path = w_total;
end
