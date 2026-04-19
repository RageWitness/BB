function GeometryCache = precompute_geometry_cache(grid_xy, ap_xy, buildings)
% PRECOMPUTE_GEOMETRY_CACHE  几何量与频率无关，预计算 (g,m) 网格
%
%   输入：
%     grid_xy    - (G x 2) 网格坐标
%     ap_xy      - (M x 2) AP 坐标
%     buildings  - 建筑列表
%
%   输出 GeometryCache：
%     d            - (G x M) 直线距离
%     Nw           - (G x M) 直线穿墙数
%     Nb           - (G x M) 阻挡建筑数
%     dc           - (G x M) 绕射路径长度（无遮挡时填 d）
%     Nc           - (G x M) 绕射拐角数
%     Nw_dif       - (G x M) 绕射沿途穿墙数
%     diff_ok      - (G x M) 是否找到绕射路径

    G = size(grid_xy, 1);
    M = size(ap_xy, 1);

    fprintf('[Geom] 预构建静态角点可见图...\n');
    t0 = tic;
    VG = build_static_corner_graph(buildings);
    fprintf('[Geom] 角点数=%d, 角点-角点边构建耗时 %.2f s\n', VG.Nc, toc(t0));

    GeometryCache.d       = zeros(G, M);
    GeometryCache.Nw      = zeros(G, M);
    GeometryCache.Nb      = zeros(G, M);
    GeometryCache.dc      = zeros(G, M);
    GeometryCache.Nc      = zeros(G, M);
    GeometryCache.Nw_dif  = zeros(G, M);
    GeometryCache.diff_ok = false(G, M);

    fprintf('[Geom] 预计算 (g,m) 几何量, G=%d, M=%d (%d 条链路)...\n', G, M, G*M);
    t1 = tic;
    print_step = max(1, floor(G / 20));

    for g = 1:G
        for m = 1:M
            p1 = grid_xy(g, :);
            p2 = ap_xy(m, :);
            d = norm(p2 - p1);
            GeometryCache.d(g, m) = d;

            [Nw, Nb] = count_wall_crossings(p1, p2, buildings);
            GeometryCache.Nw(g, m) = Nw;
            GeometryCache.Nb(g, m) = Nb;

            if Nw == 0
                GeometryCache.dc(g, m) = d;
                GeometryCache.diff_ok(g, m) = true;
            else
                pi_d = find_diffraction_path(p1, p2, buildings, VG);
                if pi_d.ok
                    GeometryCache.dc(g, m)      = pi_d.total_length;
                    GeometryCache.Nc(g, m)      = pi_d.num_corners;
                    GeometryCache.Nw_dif(g, m)  = pi_d.wall_crossings_along_path;
                    GeometryCache.diff_ok(g, m) = true;
                end
            end
        end
        if mod(g, print_step) == 0
            fprintf('  进度 %d/%d (%.1f%%), 已用 %.1f s\n', g, G, 100*g/G, toc(t1));
        end
    end
    fprintf('[Geom] 几何缓存完成, 总耗时 %.2f s\n', toc(t1));
end
