function [wall_count, blocking_building_count, hit_buildings] = count_wall_crossings(p1, p2, buildings, eps_geom)
% COUNT_WALL_CROSSINGS  线段 p1->p2 穿越所有矩形建筑边界的次数
%
%   规则：
%     - 每穿过一次矩形边界 = 1 堵墙
%     - 进入算 1 次，离开算 1 次
%     - 若 p1 在某建筑内部，离开该建筑的边界穿越要计入（自动满足，因为线段会穿出）
%     - 若 p2 在某建筑内部，进入该建筑的边界穿越要计入（同上）
%     - 仅触边不进入：不计
%     - 穿过角点：去重，不重复计数
%
%   输入：
%       p1, p2     - (1x2) 端点
%       buildings  - struct array, 含 xmin/xmax/ymin/ymax
%       eps_geom   - 浮点容差，默认 1e-9
%
%   输出：
%       wall_count               - 总边界穿越次数
%       blocking_building_count  - 实际阻挡（线段进入内部）的建筑个数
%       hit_buildings            - 1xN logical, 哪些建筑被线段穿入

    if nargin < 4 || isempty(eps_geom)
        eps_geom = 1e-9;
    end

    N = numel(buildings);
    wall_count = 0;
    blocking_building_count = 0;
    hit_buildings = false(1, N);

    p1 = p1(:)';  p2 = p2(:)';
    if norm(p2 - p1) < eps_geom
        return;
    end

    for k = 1:N
        bld = buildings(k);
        n_cross = segment_rect_boundary_crossings(p1, p2, bld, eps_geom);
        wall_count = wall_count + n_cross;
        if n_cross > 0
            % 是否真正进入内部（而非贴边）
            % 取若干内插点判断
            entered = does_segment_enter_rect(p1, p2, bld, eps_geom);
            if entered
                hit_buildings(k) = true;
                blocking_building_count = blocking_building_count + 1;
            end
        end
    end
end


function n_cross = segment_rect_boundary_crossings(p1, p2, bld, eps_geom)
% 计算线段与矩形 4 条边的"穿越"次数（去重角点）
    edges = {
        [bld.xmin bld.ymin], [bld.xmax bld.ymin];   % bottom
        [bld.xmax bld.ymin], [bld.xmax bld.ymax];   % right
        [bld.xmax bld.ymax], [bld.xmin bld.ymax];   % top
        [bld.xmin bld.ymax], [bld.xmin bld.ymin]    % left
    };

    intersect_pts = zeros(0, 2);
    for e = 1:4
        [hit, pt] = segments_intersect(p1, p2, edges{e,1}, edges{e,2}, eps_geom);
        if hit
            intersect_pts(end+1, :) = pt; %#ok<AGROW>
        end
    end

    if isempty(intersect_pts)
        n_cross = 0;
        return;
    end

    % 去重（角点会被两条边同时报）
    unique_pts = intersect_pts(1, :);
    for i = 2:size(intersect_pts, 1)
        d = sqrt(sum((unique_pts - intersect_pts(i, :)).^2, 2));
        if all(d > 1e-6)
            unique_pts(end+1, :) = intersect_pts(i, :); %#ok<AGROW>
        end
    end

    n_cross = size(unique_pts, 1);

    % 若线段实际未进入矩形内部（擦边/沿边），不计墙损
    if n_cross > 0 && ~does_segment_enter_rect(p1, p2, bld, eps_geom)
        n_cross = 0;
    end
end


function tf = does_segment_enter_rect(p1, p2, bld, eps_geom)
% 在线段上取若干内插点，检测是否有点严格位于矩形内部
    ts = linspace(0.01, 0.99, 11);
    tf = false;
    for i = 1:numel(ts)
        q = p1 + ts(i) * (p2 - p1);
        if (q(1) > bld.xmin + eps_geom) && (q(1) < bld.xmax - eps_geom) && ...
           (q(2) > bld.ymin + eps_geom) && (q(2) < bld.ymax - eps_geom)
            tf = true;
            return;
        end
    end
    % 端点在内部
    if point_in_rect(p1, bld, eps_geom) || point_in_rect(p2, bld, eps_geom)
        tf = true;
    end
end


function tf = point_in_rect(p, bld, eps_geom)
    tf = (p(1) > bld.xmin + eps_geom) && (p(1) < bld.xmax - eps_geom) && ...
         (p(2) > bld.ymin + eps_geom) && (p(2) < bld.ymax - eps_geom);
end


function [hit, pt] = segments_intersect(p1, p2, q1, q2, eps_geom)
% 二维线段相交。共线相交时返回中点；返回是否有交点和交点坐标
    hit = false;
    pt = [NaN NaN];

    r = p2 - p1;
    s = q2 - q1;
    denom = r(1)*s(2) - r(2)*s(1);
    qp = q1 - p1;

    if abs(denom) < eps_geom
        % 平行或共线
        cross_qp_r = qp(1)*r(2) - qp(2)*r(1);
        if abs(cross_qp_r) > eps_geom
            return;  % 平行不共线
        end
        % 共线：检查参数重叠
        rr = dot(r, r);
        if rr < eps_geom, return; end
        t0 = dot(qp, r) / rr;
        t1 = t0 + dot(s, r) / rr;
        tmin = max(0, min(t0, t1));
        tmax = min(1, max(t0, t1));
        if tmax + eps_geom < tmin
            return;
        end
        tmid = 0.5 * (tmin + tmax);
        pt = p1 + tmid * r;
        hit = true;
        return;
    end

    t = (qp(1)*s(2) - qp(2)*s(1)) / denom;
    u = (qp(1)*r(2) - qp(2)*r(1)) / denom;

    if t >= -eps_geom && t <= 1 + eps_geom && u >= -eps_geom && u <= 1 + eps_geom
        hit = true;
        pt = p1 + t * r;
    end
end
