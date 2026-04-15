function SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config)
% BUILD_SPATIAL_FP_SINGLE_SOURCE  构建每个频带的单源空间指纹库
%
%   SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config)
%
%   对每个频带 b、每个有效参考点 g，放置参考单源，
%   用信道模型（含 Monte Carlo 阴影采样）计算各 AP 的 RSS 指纹。
%
%   若 Config.m25.prob_shape.enable = true，还会额外计算：
%     - F_shape_prob_mu   (M x G)  — MC shape 样本均值
%     - F_shape_prob_var  (M x G)  — MC shape 样本方差
%     - ap_weight_global  (M x 1)  — 全局 AP 判别权重
%
%   输出：
%       SpatialFP.B                    - 频带数
%       SpatialFP.M                    - AP 数
%       SpatialFP.G                    - 有效网格点数
%       SpatialFP.grid_xy              - (G x 2)
%       SpatialFP.ref_power_dBm        - (1 x B) 各频带参考发射功率
%       SpatialFP.band(b).F_dBm        - (M x G)
%       SpatialFP.band(b).F_lin        - (M x G)
%       SpatialFP.band(b).std_dBm      - (M x G)
%       SpatialFP.band(b).std_lin      - (M x G)
%       SpatialFP.band(b).mean_dBm     - (1 x G)
%       SpatialFP.band(b).centered_dBm - (M x G)
%       SpatialFP.band(b).mean_lin     - (1 x G)
%       SpatialFP.band(b).centered_lin - (M x G)
%       SpatialFP.band(b).norm_l1      - (1 x G)
%       SpatialFP.band(b).F_shape_l1   - (M x G)
%       SpatialFP.band(b).F_shape_prob_mu  - (M x G) [若启用]
%       SpatialFP.band(b).F_shape_prob_var - (M x G) [若启用]
%       SpatialFP.band(b).ap_weight_global - (M x 1) [若启用]
%       SpatialFP.band(b).fc_Hz
%       SpatialFP.band(b).bw_Hz
%       SpatialFP.band(b).model

    B = Bands.B;
    M = APs.num;
    G = GridValid.Nvalid;

    % 默认参考功率
    if isfield(Config, 'm25') && isfield(Config.m25, 'ref_power_dBm')
        ref_power = Config.m25.ref_power_dBm;
    else
        ref_power = 100 * ones(1, B);  % 默认 100 dBm
    end

    % prob_shape 配置
    ps_cfg = fill_prob_shape_defaults(Config);

    SpatialFP.B = B;
    SpatialFP.M = M;
    SpatialFP.G = G;
    SpatialFP.grid_xy       = GridValid.xy;
    SpatialFP.ref_power_dBm = ref_power;

    fprintf('[M2.5] 构建 SpatialFP: M=%d AP, G=%d 网格点, B=%d 频带\n', M, G, B);
    if ps_cfg.enable
        fprintf('[M2.5] 概率 shape 指纹已启用\n');
    end

    for b = 1:B
        F_dBm   = zeros(M, G);
        F_lin   = zeros(M, G);
        Std_dBm = zeros(M, G);
        Std_lin = zeros(M, G);

        % 若启用 prob_shape，需要收集 MC 样本级 shape 统计
        if ps_cfg.enable
            shape_mu_acc  = zeros(M, G);   % 累加 shape 样本均值
            shape_var_acc = zeros(M, G);   % 累加 shape 样本方差
        end

        for g = 1:G
            if ps_cfg.enable
                % 请求第 5 个输出：MC 线性功率样本矩阵
                [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin, rss_lin_samp] = ...
                    simulate_reference_fp_point_m1bridge( ...
                        GridValid.xy(g, :), b, ref_power(b), APs, Bands, Config);

                % 对每个 MC 样本做 L1 归一化得到 shape 样本
                shape_samples = rss_lin_samp ./ (sum(rss_lin_samp, 1) + ps_cfg.eps_shape);  % M x N_mc

                % 逐维统计
                shape_mu_acc(:, g)  = mean(shape_samples, 2);    % M x 1
                shape_var_acc(:, g) = var(shape_samples, 0, 2);  % M x 1
            else
                [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin] = ...
                    simulate_reference_fp_point_m1bridge( ...
                        GridValid.xy(g, :), b, ref_power(b), APs, Bands, Config);
            end

            F_dBm(:, g)   = rss_dBm;
            F_lin(:, g)   = rss_lin;
            Std_dBm(:, g) = rss_std_dBm;
            Std_lin(:, g) = rss_std_lin;
        end

        % 组装单频带结构
        band_fp.F_dBm   = F_dBm;
        band_fp.F_lin   = F_lin;
        band_fp.std_dBm = Std_dBm;
        band_fp.std_lin = Std_lin;
        band_fp.fc_Hz   = Bands.fc_Hz(b);
        band_fp.bw_Hz   = Bands.bw_Hz(b);
        band_fp.model   = Bands.model{b};

        % 预计算均值和中心化指纹（dBm + 线性 两域）+ L1 shape
        band_fp = precompute_centered_fp_wknn(band_fp);

        % ---- 概率 shape 指纹 & 全局 AP 权重 ----
        if ps_cfg.enable
            band_fp.F_shape_prob_mu  = shape_mu_acc;    % M x G
            band_fp.F_shape_prob_var = shape_var_acc;    % M x G

            % 全局 AP 判别权重：
            %   w_m = Var_g(log(phi_{g,m} + eps)) / (Mean_g(sigma_{g,m}^2) + eps)
            log_phi = log(F_lin + ps_cfg.eps_weight);                  % M x G
            var_across_grid    = var(log_phi, 0, 2);                   % M x 1
            mean_var_within    = mean(Std_lin.^2, 2);                  % M x 1
            w_raw = var_across_grid ./ (mean_var_within + ps_cfg.eps_weight);  % M x 1

            % 归一化
            switch ps_cfg.weight_normalize
                case 'mean1'
                    w_norm = w_raw / (mean(w_raw) + ps_cfg.eps_weight);
                case 'none'
                    w_norm = w_raw;
                otherwise
                    w_norm = w_raw / (mean(w_raw) + ps_cfg.eps_weight);
            end

            % 裁剪
            if ps_cfg.weight_clip_enable
                w_norm = min(max(w_norm, ps_cfg.weight_clip_min), ps_cfg.weight_clip_max);
            end

            band_fp.ap_weight_global = w_norm;  % M x 1

            fprintf('  Band %d: AP 权重 [%.3f, %.3f], shape_var [%.2e, %.2e]\n', ...
                b, min(w_norm), max(w_norm), ...
                min(shape_var_acc(:)), max(shape_var_acc(:)));
        end

        SpatialFP.band(b) = band_fp;

        fprintf('  Band %d (%s): F_dBm [%.1f, %.1f], F_lin [%.2e, %.2e], shape_l1 [%.4f, %.4f]\n', ...
            b, Bands.name{b}, min(F_dBm(:)), max(F_dBm(:)), ...
            min(F_lin(:)), max(F_lin(:)), ...
            min(band_fp.F_shape_l1(:)), max(band_fp.F_shape_l1(:)));
    end

    % ---- 离线 hard-negative 构建（需要概率型 shape 已就绪） ----
    hn_cfg = fill_hardneg_defaults(Config);
    if hn_cfg.enable
        fprintf('[M2.5] 构建 hard-negative 集合...\n');
        for b = 1:B
            hn_out = build_hardneg_for_band( ...
                SpatialFP.band(b), SpatialFP.grid_xy, hn_cfg);
            % 逐字段赋值，避免 struct array 字段不一致
            SpatialFP.band(b).hardneg_idx        = hn_out.hardneg_idx;
            SpatialFP.band(b).hardneg_shape_dist = hn_out.hardneg_shape_dist;
            SpatialFP.band(b).hardneg_phys_dist  = hn_out.hardneg_phys_dist;
            SpatialFP.band(b).ap_weight_hn       = hn_out.ap_weight_hn;
            fprintf('  Band %d: hardneg 完成, ap_weight_hn [%.3f, %.3f]\n', ...
                b, min(SpatialFP.band(b).ap_weight_hn(:)), ...
                max(SpatialFP.band(b).ap_weight_hn(:)));
        end
    end

    fprintf('[M2.5] SpatialFP 构建完成\n');
end


%% ==================== 局部函数 ====================

function ps_cfg = fill_prob_shape_defaults(Config)
% FILL_PROB_SHAPE_DEFAULTS  填充 M2.5 prob_shape 默认配置

    ps_cfg = struct();
    ps_cfg.enable            = true;
    ps_cfg.eps_shape         = 1e-12;
    ps_cfg.eps_weight        = 1e-12;
    ps_cfg.weight_clip_enable = true;
    ps_cfg.weight_clip_min   = 0.25;
    ps_cfg.weight_clip_max   = 4.0;
    ps_cfg.weight_normalize  = 'mean1';

    if isfield(Config, 'm25') && isfield(Config.m25, 'prob_shape')
        ps_in = Config.m25.prob_shape;
        fns = fieldnames(ps_cfg);
        for i = 1:numel(fns)
            if isfield(ps_in, fns{i})
                ps_cfg.(fns{i}) = ps_in.(fns{i});
            end
        end
    end
end


function hn_cfg = fill_hardneg_defaults(Config)
% FILL_HARDNEG_DEFAULTS  填充 M2.5 hardneg 默认配置

    hn_cfg = struct();
    hn_cfg.enable            = true;
    hn_cfg.source_shape      = 'prob_mu';   % 'prob_mu' or 'shape_l1'
    hn_cfg.r_far_m           = 25;
    hn_cfg.top_conf          = 6;
    hn_cfg.eps               = 1e-12;
    hn_cfg.weight_clip_enable = true;
    hn_cfg.weight_clip_min   = 0.25;
    hn_cfg.weight_clip_max   = 4.0;
    hn_cfg.k_probe           = 2500;         % shape KNN 初始探测数（需 > π·r_far²/grid_step² ≈ 1963）
    hn_cfg.k_probe_max       = 5000;         % 自适应扩展上限
    hn_cfg.block_size        = 2000;        % 分块大小

    if isfield(Config, 'm25') && isfield(Config.m25, 'hardneg')
        hn_in = Config.m25.hardneg;
        fns = fieldnames(hn_cfg);
        for i = 1:numel(fns)
            if isfield(hn_in, fns{i})
                hn_cfg.(fns{i}) = hn_in.(fns{i});
            end
        end
    end
end


function band_fp = build_hardneg_for_band(band_fp, grid_xy, hn_cfg)
% BUILD_HARDNEG_FOR_BAND  稀疏 hard-negative 挖掘（无 G×G 矩阵）
%
%   对每个网格点 g：
%     1. 在 shape 空间中分块查询前 K_probe 个最近邻候选
%     2. 对候选计算物理距离，筛选 d_phys >= r_far_m
%     3. 若满足条件的候选不足 top_conf，自适应扩大至 K_probe_max
%     4. 按 shape 距离从小到大保留前 top_conf 个作为 hard-negative
%     5. 计算局部 AP 判别权重 w_hn(:,g)
%
%   严禁创建 D_shape(G,G) 或 D_phys(G,G)。
%   峰值内存 ~ block_size * G * 8 bytes（shape 距离分块）。

    [M, G] = size(band_fp.F_lin);

    % ---- 选择 shape 源 ----
    if strcmp(hn_cfg.source_shape, 'prob_mu') && isfield(band_fp, 'F_shape_prob_mu')
        S = band_fp.F_shape_prob_mu;    % M x G
    else
        S = band_fp.F_shape_l1;         % M x G  回退
    end

    % 获取 shape 方差（若有）
    if isfield(band_fp, 'F_shape_prob_var')
        V = band_fp.F_shape_prob_var;    % M x G
    else
        V = zeros(M, G);
    end

    % ---- 预计算 ||S(:,g)||^2 ----
    S_sq_all = sum(S.^2, 1);            % 1 x G

    top_conf   = hn_cfg.top_conf;
    r_far      = hn_cfg.r_far_m;
    eps_hn     = hn_cfg.eps;
    K_probe0   = min(hn_cfg.k_probe, G - 1);
    K_probe_mx = min(hn_cfg.k_probe_max, G - 1);
    BLK        = min(hn_cfg.block_size, G);

    % ---- 稀疏结果 ----
    hardneg_idx        = zeros(G, top_conf);
    hardneg_shape_dist = zeros(G, top_conf);
    hardneg_phys_dist  = zeros(G, top_conf);
    ap_weight_hn       = ones(M, G);

    xy_all = grid_xy;                    % G x 2

    % ================================================================
    %  第一遍：分块计算 shape 距离，取每个 g 的 top-K_probe0 近邻索引
    %  存储：knn_idx (G x K_probe0) + knn_dist (G x K_probe0)
    %  峰值内存 ~ BLK * G * 8 bytes
    % ================================================================
    knn_idx  = zeros(G, K_probe0);
    knn_dist = zeros(G, K_probe0);

    for blk_s = 1:BLK:G
        blk_e = min(blk_s + BLK - 1, G);
        idx_b = blk_s:blk_e;
        Nb = numel(idx_b);

        % shape L2 距离 Nb x G（分块展开，不缓存 G×G）
        S_blk = S(:, idx_b);                          % M x Nb
        D2 = S_sq_all(idx_b)' + S_sq_all ...
             - 2 * (S_blk' * S);                       % Nb x G
        D2 = max(D2, 0);

        % 排除自身：设为 inf
        for i = 1:Nb
            D2(i, idx_b(i)) = inf;
        end

        % partial sort：取 top K_probe0 最小
        [sorted_d2, sorted_j] = sort(D2, 2, 'ascend'); % Nb x G  → 取前 K 列
        knn_idx(idx_b, :)  = sorted_j(:, 1:K_probe0);
        knn_dist(idx_b, :) = sqrt(sorted_d2(:, 1:K_probe0));
    end

    % ================================================================
    %  第二遍：逐点筛选 + 自适应扩展
    % ================================================================
    n_expand = 0;  % 统计自适应扩展次数

    for g = 1:G
        [hn_g, d_shape_g, d_phys_g, need_expand] = select_hardneg_from_knn( ...
            g, knn_idx(g,:), knn_dist(g,:), xy_all, r_far, top_conf);

        % ---- 自适应扩展 ----
        if need_expand && K_probe_mx > K_probe0
            % 只对当前点 g 计算全库 shape 距离（1 x G 向量）
            d2_g = S_sq_all(g) + S_sq_all - 2 * (S(:,g)' * S);  % 1 x G
            d2_g = max(d2_g, 0);
            d2_g(g) = inf;

            [sorted_d2_g, sorted_j_g] = sort(d2_g, 'ascend');
            K_ext = min(K_probe_mx, numel(sorted_j_g));
            ext_idx  = sorted_j_g(1:K_ext);
            ext_dist = sqrt(sorted_d2_g(1:K_ext));

            [hn_g, d_shape_g, d_phys_g, ~] = select_hardneg_from_knn( ...
                g, ext_idx, ext_dist, xy_all, r_far, top_conf);
            n_expand = n_expand + 1;
        end

        n_take = numel(hn_g);
        if n_take > 0
            hardneg_idx(g, 1:n_take)        = hn_g;
            hardneg_shape_dist(g, 1:n_take) = d_shape_g;
            hardneg_phys_dist(g, 1:n_take)  = d_phys_g;

            % ---- 局部 AP 判别权重 ----
            mu_g  = S(:, g);
            var_g = V(:, g);
            mu_H  = S(:, hn_g);
            var_H = V(:, hn_g);

            diff_abs = abs(mu_g - mu_H);
            denom    = sqrt(var_g + var_H + eps_hn);
            a_gm     = mean(diff_abs ./ denom, 2);       % M x 1

            w_g = a_gm / (mean(a_gm) + eps_hn);
            if hn_cfg.weight_clip_enable
                w_g = min(max(w_g, hn_cfg.weight_clip_min), hn_cfg.weight_clip_max);
            end
            ap_weight_hn(:, g) = w_g;
        end
    end

    if n_expand > 0
        fprintf('    自适应扩展 K_probe: %d / %d 点\n', n_expand, G);
    end

    band_fp.hardneg_idx        = hardneg_idx;         % G x top_conf
    band_fp.hardneg_shape_dist = hardneg_shape_dist;   % G x top_conf
    band_fp.hardneg_phys_dist  = hardneg_phys_dist;    % G x top_conf
    band_fp.ap_weight_hn       = ap_weight_hn;          % M x G
end


function [hn_idx, d_shape, d_phys, need_expand] = select_hardneg_from_knn( ...
        g, cand_idx, cand_dist, xy_all, r_far, top_conf)
% SELECT_HARDNEG_FROM_KNN  从 KNN 候选中筛选满足物理距离约束的 hard-neg
%
%   输入：
%       g          - 查询点全局索引
%       cand_idx   - (1 x K) KNN 候选全局索引
%       cand_dist  - (1 x K) 对应 shape L2 距离
%       xy_all     - (G x 2) 全库坐标
%       r_far      - 物理距离阈值
%       top_conf   - 需要保留的 hard-neg 数量
%
%   输出：
%       hn_idx      - 满足条件的 hard-neg 全局索引（最多 top_conf 个）
%       d_shape     - 对应 shape 距离
%       d_phys      - 对应物理距离
%       need_expand - 是否候选不足需要扩展

    % 计算候选的物理距离（仅 K_probe 个点，非 G 个）
    dx = xy_all(g, 1) - xy_all(cand_idx, 1);
    dy = xy_all(g, 2) - xy_all(cand_idx, 2);
    phys_d = sqrt(dx.^2 + dy.^2);                    % K x 1

    % 筛选物理距离 >= r_far
    far_mask = phys_d >= r_far;
    far_pos = find(far_mask);

    if isempty(far_pos)
        hn_idx = []; d_shape = []; d_phys = [];
        need_expand = true;
        return;
    end

    % cand_dist 已按 shape 距离升序排列，far_pos 中最小的就是最近
    n_take = min(top_conf, numel(far_pos));
    % far_pos 中的元素对应 cand_dist 的位置，取 shape 距离最小的 n_take 个
    [~, re_order] = sort(cand_dist(far_pos), 'ascend');
    sel = far_pos(re_order(1:n_take));

    hn_idx  = cand_idx(sel);
    d_shape = cand_dist(sel);
    d_phys  = phys_d(sel);

    need_expand = (n_take < top_conf);
end
