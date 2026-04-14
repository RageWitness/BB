function scr = screen_candidates_by_area_m4( ...
    d_shape_vec, grid_xy, scr_cfg)
% SCREEN_CANDIDATES_BY_AREA_M4  物理空间区域筛选（Area Screening）
%
%   scr = screen_candidates_by_area_m4(d_shape_vec, grid_xy, scr_cfg)
%
%   d_shape 只用于排序选出预候选；之后的验证圆完全定义在物理坐标空间。
%
%   算法：
%     1. 按 d_shape 选出 K_shape 个预候选（shape 空间排序）
%     ---- 以下全部在物理坐标空间 ----
%     2. 取前两名 shape 邻居的物理坐标 x_g1, x_g2
%     3. 物理距离 d12_phys = ||x_g1 - x_g2||_2
%        若 d12_phys > d_pair_max → 回退
%     4. 验证圆心 c = (x_g1 + x_g2) / 2   （物理中点）
%     5. 自适应半径 r = min(max(r_min, alpha * d12_phys), r_max)  （物理距离）
%     6. 对每个预候选 g，计算其物理坐标到圆心的欧氏距离
%        保留 ||x_g - c||_2 <= r 的候选
%     7. 若通过筛选的候选 < k_min_valid → 回退
%
%   输入：
%       d_shape_vec - (1 x G) 全库形状距离（仅用于排序选预候选）
%       grid_xy     - (G x 2) 网格点物理坐标 (m)
%       scr_cfg     - 配置结构体
%
%   输出：
%       scr - 筛选结果结构体（所有距离/坐标均为物理空间量）

    G = numel(d_shape_vec);

    % ---- 初始化输出（回退默认值） ----
    scr = struct();
    scr.screen_applied     = false;
    scr.fallback_used      = false;
    scr.fallback_reason    = '';
    scr.pre_idx_shape      = [];        % shape 排序得到的预候选全局索引
    scr.pre_dist_shape     = [];        % 对应的 d_shape 值（仅记录，不参与后续）
    scr.pre_xy_phys        = [];        % 预候选的物理坐标 (K_shape x 2)
    scr.center_xy_phys     = [NaN, NaN];% 验证圆心（物理坐标中点）
    scr.radius_phys_m      = NaN;       % 验证半径（物理距离）
    scr.pair_dist12_phys_m = NaN;       % 前两名的物理间距
    scr.dist_to_center_phys_m = [];     % 每个预候选到圆心的物理距离
    scr.valid_mask_phys    = [];        % 物理圆内/外判定 (1 x K_shape logical)
    scr.valid_idx          = [];        % 通过筛选的全局网格索引
    scr.rejected_idx       = [];        % 被剔除的全局网格索引
    scr.valid_count        = 0;

    % ---- 检查是否启用 ----
    if ~scr_cfg.enable
        scr.fallback_used   = true;
        scr.fallback_reason = 'disabled';
        return;
    end

    %% ===== Step 1: 按 d_shape 排序选预候选 =====
    % （d_shape 的唯一用途：排序。之后不再参与任何圆的计算。）
    K_shape = min(scr_cfg.k_shape, G);
    [sorted_dist, sorted_idx] = sort(d_shape_vec, 'ascend');

    pre_idx_shape  = sorted_idx(1:K_shape);     % 预候选全局索引
    pre_dist_shape = sorted_dist(1:K_shape);     % 对应 d_shape（仅记录）
    pre_xy_phys    = grid_xy(pre_idx_shape, :);  % 预候选物理坐标 (K_shape x 2)

    scr.pre_idx_shape  = pre_idx_shape;
    scr.pre_dist_shape = pre_dist_shape;
    scr.pre_xy_phys    = pre_xy_phys;

    % ---- 保护: 预候选不足 ----
    if K_shape < 2
        scr.fallback_used   = true;
        scr.fallback_reason = 'pre_candidates_less_than_2';
        return;
    end

    %% ===== Step 2: 前两名 shape 邻居的物理坐标 =====
    x_g1 = pre_xy_phys(1, :);   % 1 x 2，物理坐标
    x_g2 = pre_xy_phys(2, :);   % 1 x 2，物理坐标

    %% ===== Step 3: 前两名的物理距离 =====
    pair_dist12_phys_m = norm(x_g1 - x_g2, 2);
    scr.pair_dist12_phys_m = pair_dist12_phys_m;

    % ---- 保护: 前两名物理距离过远 ----
    if pair_dist12_phys_m > scr_cfg.d_pair_max_m
        scr.fallback_used   = true;
        scr.fallback_reason = sprintf( ...
            'pair_phys_dist_%.1fm_exceeds_max_%.1fm', ...
            pair_dist12_phys_m, scr_cfg.d_pair_max_m);
        return;
    end

    %% ===== Step 4: 验证圆心 = 物理坐标中点 =====
    center_xy_phys = (x_g1 + x_g2) / 2;
    scr.center_xy_phys = center_xy_phys;

    %% ===== Step 5: 自适应半径（基于物理距离） =====
    %   r = min( max(r_min, alpha * d12_phys), r_max )
    radius_phys_m = min( ...
        max(scr_cfg.r_min_m, scr_cfg.alpha * pair_dist12_phys_m), ...
        scr_cfg.r_max_m );
    scr.radius_phys_m = radius_phys_m;

    %% ===== Step 6: 物理空间圆内判定 =====
    % 每个预候选到圆心的物理欧氏距离
    dist_to_center_phys_m = sqrt(sum( ...
        (pre_xy_phys - center_xy_phys).^2, 2))';   % 1 x K_shape
    scr.dist_to_center_phys_m = dist_to_center_phys_m;

    % 保留条件：物理距离 <= 验证半径
    valid_mask_phys = dist_to_center_phys_m <= radius_phys_m;  % 1 x K_shape
    scr.valid_mask_phys = valid_mask_phys;
    scr.valid_idx    = pre_idx_shape(valid_mask_phys);
    scr.rejected_idx = pre_idx_shape(~valid_mask_phys);
    scr.valid_count  = sum(valid_mask_phys);

    % ---- 保护: 筛选后剩余太少 ----
    if scr.valid_count < scr_cfg.k_min_valid
        scr.fallback_used   = true;
        scr.fallback_reason = sprintf( ...
            'valid_count_%d_below_min_%d', ...
            scr.valid_count, scr_cfg.k_min_valid);
        return;
    end

    %% ===== 筛选成功 =====
    scr.screen_applied = true;
end
