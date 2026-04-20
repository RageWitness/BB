%% MAIN_SIMULATION  主仿真入口 — 新框架
%
%  流程：
%    1. M0 初始化 + M1 初始化 + M2.5 指纹库构建
%    2. 主循环：逐帧调度源 -> 生成 RSS 观测
%    3. 外部源类型标签输入 + 帧级联动
%    4. M4 WKNN 定位（基于外部标签而非 M3 分类）
%    5. 硬约束验证
%    6. 统计输出 + 可视化
%
%  注意：M3 分类模块已从主链路断开。
%        源类型由外部输入（当前从 M0 真值日志提取）。

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  频谱指纹匹配机遇式增强定位 — 主仿真（新框架）\n');
fprintf('============================================\n\n');

%% ========== 第一阶段：初始化 ==========

% --- M0 初始化 ---
sim_override.ap.layout = 'explicit';
sim_override.ap.pos_xy = [
   65.2850   ,13.7581;
   14.0990 , 186.6792;
  192.2744  ,255.9371;
  203.5478  , 25.5561;
  273.7034  , 98.0424;
   67.5742,  113.4034;
  288.0000  ,215.4216;
   75.6478  ,274.5550
];
sim_override.fp.grid_step = 3;   % 3m grid

% --- 指纹库缓存开关 ---
fp_cache_enable = true;                          % true: 命中即加载、未命中则构建并保存
fp_cache_file   = 'cache/SpatialFP_mwm.mat';

% --- M4 匹配指纹类型选择 ---
sim_override.m4.fingerprint_type = 'centered_dBm';   % rf_minmax | rf_raw | shape_l1 | centered_dBm | legacy
sim_override.m4.fp_distance      = 'L2';          % L1 | L2

% sim_override.m4.distance_mode = 'shape_scale';

sim_override.debug.expose_true_source_state = true;


[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);

% --- M1 初始化 ---
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

% --- M2.5 指纹库（新框架：主模式 rf_minmax + legacy 兼容）---
if fp_cache_enable && exist(fp_cache_file, 'file')
    fprintf('[Cache] 加载指纹库缓存: %s\n', fp_cache_file);
    S = load(fp_cache_file);
    SpatialFP    = S.SpatialFP;
    SignatureLib = S.SignatureLib;
else
    [SpatialFP, SignatureLib] = init_m25_single_source_fp( ...
        APs, Bands, GridValid, Config, SourceTemplates);
    if fp_cache_enable
        cache_dir = fileparts(fp_cache_file);
        if ~isempty(cache_dir) && ~exist(cache_dir, 'dir'), mkdir(cache_dir); end
        save(fp_cache_file, 'SpatialFP', 'SignatureLib', '-v7.3');
        fprintf('[Cache] 指纹库已保存: %s\n', fp_cache_file);
    end
end

% --- 全局参数 ---
T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;
dt = Config.m0.dt;

fprintf('\n--- 仿真参数 ---\n');
fprintf('  区域: [%.0f, %.0f] x [%.0f, %.0f] m\n', Config.area.x_range, Config.area.y_range);
fprintf('  AP: %d 个 | 有效网格: %d 点\n', M, GridValid.Nvalid);
fprintf('  频带: %d 个 | 总帧数: %d | dt=%.1f s\n', B, T, dt);
fprintf('  指纹模式: %s\n', SpatialFP.fingerprint_mode);
fprintf('  M4 指纹匹配: type=%s, dist=%s\n', Config.m4.fingerprint_type, Config.m4.fp_distance);
for b = 1:B
    fprintf('  Band %d: %s, fc=%.1f MHz, model=%s\n', ...
        b, Bands.name{b}, Bands.fc_Hz(b)/1e6, Bands.model{b});
end

%% ========== 第二阶段：主循环 ==========

fprintf('\n--- 主循环开始 (%d 帧) ---\n', T);

FrameStates      = cell(1, T);
ObsFrames        = cell(1, T);
Y_dBm_all        = zeros(M, B, T);
Y_lin_all        = zeros(M, B, T);
occupancy_matrix = zeros(T, B);

tic;
for t = 1:T
    % --- M0：源调度 ---
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);

    % --- M1：信道观测 ---
    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);

    % --- 存储 ---
    FrameStates{t} = FrameState_t;
    ObsFrames{t}   = ObsFrame;
    Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;
    Y_lin_all(:, :, t)  = ObsFrame.Y_lin;

    % 记录占用类型
    for b = 1:B
        apb = FrameState_t.active_per_band(b);
        if apb.has_source
            switch apb.source_type
                case 'broadband_cal',  occupancy_matrix(t, b) = 1;
                case 'opportunistic',  occupancy_matrix(t, b) = 2;
                case 'target',         occupancy_matrix(t, b) = 4;
                case 'persistent_cal', occupancy_matrix(t, b) = 3;
            end
        end
    end
end
elapsed = toc;
fprintf('主循环完成: %d 帧, 耗时 %.2f s (%.1f 帧/s)\n', T, elapsed, T/elapsed);

%% ========== 第 2.5a 阶段：构建 DrivenSourceInput ==========

fprintf('\n--- 构建 DrivenSourceInput ---\n');
DSI_all = cell(T, B);
for t = 1:T
    for b = 1:B
        DSI_all{t, b} = build_driven_source_input( ...
            FrameStates{t}.active_per_band(b), Config);
    end
end
n_dsi_active = sum(cellfun(@(d) d.label > 0, DSI_all), 'all');
fprintf('[DSI] 构建完成: %d x %d 格, 有源 %d 格\n', T, B, n_dsi_active);

%% ========== 第 2.5b 阶段：外部源类型输入 + M4 定位 ==========
%
%  新框架：M3 分类已断开，源类型由外部输入。
%  当前从 M0 真值日志提取外部标签（仿真环境下的 oracle 模式）。
%  后续对接真实系统时，改为从外部接口读取。

fprintf('\n--- 外部源类型标签输入 ---\n');

% --- 从真值日志构建外部标签 ---
ExternalLabelsRaw = build_labels_from_truth_log(M0Logs, Config);
LabelTable = ingest_external_source_labels(ExternalLabelsRaw, Config);

% --- 标签绑定到帧 ---
SourceContext = bind_external_labels_to_events(LabelTable, FrameStates, Config);

% --- M4：WKNN 定位（基于外部标签，不再依赖 M3 EventList）---
%     过渡态：M4 仍使用旧的 EventList 接口，
%     这里从外部标签构建一个兼容的 EventList 传入 M4。
EventList = build_event_list_from_context(SourceContext, Y_dBm_all, Y_lin_all, Config);

fprintf('\n');
LocResults = run_m4_wknn_localization( ...
    EventList, SpatialFP, FrameStates, Config);

%% ========== 第三阶段：验证 ==========

fprintf('\n--- 硬约束验证 ---\n');
for t = 1:T
    for b = 1:B
        n_active = double(FrameStates{t}.active_per_band(b).has_source);
        assert(n_active <= 1, 'N_active(%d,%d) = %d > 1 违反硬约束！', b, t, n_active);
    end
end
fprintf('[PASS] 所有频带所有帧 N_active(b,t) <= 1\n');

%% ========== 第四阶段：统计 ==========

fprintf('\n--- 真值日志统计 ---\n');
fprintf('  TruthLogAll   : %d 条\n', numel(M0Logs.TruthLogAll));
fprintf('  TruthLogTarget: %d 条\n', numel(M0Logs.TruthLogTarget));

% 频带占用统计
type_names = {'空闲', 'broadband_cal', 'opportunistic', 'persistent_cal', 'target'};
fprintf('\n--- 各频带占用统计 (帧数 / %d) ---\n', T);
fprintf('  %-6s', 'Band');
for ti = 1:5
    fprintf('%-18s', type_names{ti});
end
fprintf('\n');
for b = 1:B
    fprintf('  %-6d', b);
    for ti = 0:4
        cnt = sum(occupancy_matrix(:, b) == ti);
        fprintf('%-18d', cnt);
    end
    fprintf('\n');
end

% 外部标签统计
fprintf('\n--- 外部标签统计 ---\n');
fprintf('  标签总数: %d\n', SourceContext.n_labels);
if ~isempty(LabelTable)
    LC_stat = source_label_constants();
    labels_arr = [LabelTable.label];
    for lb = LC_stat.ALL_LABELS
        cnt = sum(labels_arr == lb);
        fprintf('  label=%d (%s): %d\n', lb, LC_stat.name_map(lb), cnt);
    end
end

% RSS 统计
fprintf('\n--- RSS 统计 (有源帧均值) ---\n');
for b = 1:B
    src_frames = find(occupancy_matrix(:, b) > 0);
    if ~isempty(src_frames)
        rss_src = squeeze(Y_dBm_all(:, b, src_frames));
        fprintf('  Band %d (%s): mean=%.1f dBm, std=%.1f dB, %d 帧有源\n', ...
            b, Bands.name{b}, mean(rss_src(:)), std(rss_src(:)), numel(src_frames));
    else
        fprintf('  Band %d (%s): 无有源帧\n', b, Bands.name{b});
    end
end

%% ========== 第五阶段：可视化 ==========

% ===== 图1：频带占用时间线 =====
figure('Name', '频带占用时间线', 'Position', [50 50 1200 350]);
cmap = [1 1 1;         % 0: 空闲
        0.2 0.6 1;     % 1: trusted
        0 0.8 0.4;     % 2: prior_pos
        1 0.7 0;       % 3: prior_time
        0.9 0.2 0.2];  % 4: target
imagesc(1:T, 1:B, occupancy_matrix');
colormap(cmap);
cbar = colorbar;
cbar.Ticks = [0, 1, 2, 3, 4];
cbar.TickLabels = type_names;
xlabel('帧号'); ylabel('频带');
title('M0 频带占用时间线');
set(gca, 'YTick', 1:B);
grid on;

% ===== 图2：RSS 热力图 (AP x 帧) =====
figure('Name', 'RSS 热力图', 'Position', [50 450 1200 700]);
for b = 1:B
    subplot(2, 2, b);
    imagesc(1:T, 1:M, squeeze(Y_dBm_all(:, b, :)));
    colorbar;
    xlabel('帧号'); ylabel('AP');
    title(sprintf('Band %d (%s) RSS [dBm]', b, Bands.name{b}));
    set(gca, 'YTick', 1:4:M);
end
sgtitle('M1 RSS 热力图 (AP x 帧)');

% ===== 图3：场景布局 — 各类源累积位置 =====
figure('Name', '场景布局', 'Position', [700 50 700 600]);
plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.85 0.85 0.85], 'MarkerSize', 3);
hold on;

% --- 绘制建筑（半透明灰色色块 + 深色边框）---
if isfield(Config, 'm1') && isfield(Config.m1, 'channel') && ...
        isfield(Config.m1.channel, 'buildings') && ~isempty(Config.m1.channel.buildings)
    blds = Config.m1.channel.buildings;
    for k = 1:numel(blds)
        b = blds(k);
        xv = [b.xmin b.xmax b.xmax b.xmin];
        yv = [b.ymin b.ymin b.ymax b.ymax];
        patch(xv, yv, [0.35 0.45 0.55], 'FaceAlpha', 0.25, ...
            'EdgeColor', [0.25 0.30 0.40], 'LineWidth', 0.8);
    end
end

plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

if ~isempty(M0Logs.TruthLogAll)
    all_types = {M0Logs.TruthLogAll.source_type};
    all_xy    = reshape([M0Logs.TruthLogAll.true_pos_xy], 2, [])';
    all_ids   = [M0Logs.TruthLogAll.instance_id];
    all_keys  = {M0Logs.TruthLogAll.template_key};

    mask_tr = strcmp(all_types, 'broadband_cal');
    h_tr = plot_unique_sources(all_xy, all_ids, mask_tr, 'p', [0.2 0.5 1], 16);
    mask_pp = strcmp(all_types, 'persistent_cal');
    h_pp = plot_unique_sources(all_xy, all_ids, mask_pp, 's', [0 0.8 0.4], 10);
    mask_pt = strcmp(all_types, 'opportunistic');
    h_pt = plot_unique_sources(all_xy, all_ids, mask_pt, 'd', [1 0.7 0], 9);
    mask_tg = strcmp(all_types, 'target');
    h_tg = plot_unique_sources(all_xy, all_ids, mask_tg, 'o', [0.1 0.2 0.2], 7);
end

xlabel('X (m)'); ylabel('Y (m)');
title('场景布局 — 各类源累积出现位置');
axis equal; grid on;
xlim(Config.area.x_range); ylim(Config.area.y_range);
leg_h = {}; leg_s = {};
leg_h{end+1} = plot(NaN, NaN, '.', 'Color', [0.85 0.85 0.85], 'MarkerSize', 10);
leg_s{end+1} = '有效网格点';
leg_h{end+1} = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
leg_s{end+1} = 'AP';
if ~isempty(h_tr), leg_h{end+1} = h_tr; leg_s{end+1} = 'broadband\_cal'; end
if ~isempty(h_pp), leg_h{end+1} = h_pp; leg_s{end+1} = 'persistent\_cal'; end
if ~isempty(h_pt), leg_h{end+1} = h_pt; leg_s{end+1} = 'opportunistic'; end
if ~isempty(h_tg), leg_h{end+1} = h_tg; leg_s{end+1} = 'target'; end
legend([leg_h{:}], leg_s, 'Location', 'best');

% ===== 图4：最后有源帧 RSS vs 距离 =====
figure('Name', 'RSS vs 距离', 'Position', [700 700 1000 400]);
for b = 1:B
    src_frames = find(occupancy_matrix(:, b) > 0);
    if isempty(src_frames), continue; end
    t_last = src_frames(end);
    src_xy = FrameStates{t_last}.active_per_band(b).true_pos_xy;
    dist_m = sqrt(sum((APs.pos_xy - src_xy).^2, 2));
    subplot(1, B, b);
    scatter(dist_m, Y_dBm_all(:, b, t_last), 40, 'filled');
    xlabel('距离 (m)'); ylabel('RSS (dBm)');
    title(sprintf('Band %d 帧%d', b, t_last));
    grid on;
end
sgtitle('RSS vs 源-AP 距离');


% ===== 图5：M4 定位结果散点图 =====
if ~isempty(LocResults)
    figure('Name', 'M4 定位结果', 'Position', [100 100 800 700]);
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    for k = 1:numel(LocResults)
        lr = LocResults(k);
        if isempty(lr.true_pos_xy) || all(lr.true_pos_xy == 0)
            continue;
        end
        switch lr.type_hat
            case 'target'
                clr = [0.9 0.2 0.2];
                short_name = 'target';
            case 'opportunistic'
                clr = [1.0 0.7 0.0];
                short_name = 'opp';
            otherwise
                clr = [0.5 0.5 0.5];
                short_name = lr.type_hat;
        end
        plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', ...
            'Color', clr, 'MarkerSize', 8, 'MarkerFaceColor', clr);
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
            'Color', clr*0.6, 'MarkerSize', 10, 'LineWidth', 2);
        plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
             [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
             '-', 'Color', min(1, clr*0.65 + 0.25), 'LineWidth', 1);
        text(lr.true_pos_xy(1)+2, lr.true_pos_xy(2)+2, ...
            sprintf('B%d %s', lr.band_id, short_name), ...
            'Color', clr*0.75, 'FontSize', 8, 'FontWeight', 'bold');
        text(lr.est_pos_xy(1)+2, lr.est_pos_xy(2)-4, ...
            sprintf('B%d est', lr.band_id), ...
            'Color', clr*0.55, 'FontSize', 8);
    end
    xlabel('X (m)'); ylabel('Y (m)');
    title('M4 WKNN 定位结果 (o=真值, x=估计)');
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);
    leg_h = gobjects(5, 1);
    leg_h(1) = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    leg_h(2) = plot(NaN, NaN, 'o', 'Color', [0.9 0.2 0.2], 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.2 0.2]);
    leg_h(3) = plot(NaN, NaN, 'x', 'Color', [0.9 0.2 0.2]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    leg_h(4) = plot(NaN, NaN, 'o', 'Color', [1 0.7 0], 'MarkerSize', 8, 'MarkerFaceColor', [1 0.7 0]);
    leg_h(5) = plot(NaN, NaN, 'x', 'Color', [1 0.7 0]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    legend(leg_h, {'AP', 'target true', 'target est', 'opportunistic true', 'opportunistic est'}, ...
        'Location', 'best');
end

% ===== 图6：分频带定位误差统计 =====
if ~isempty(LocResults)
    valid_mask = ~arrayfun(@(r) isnan(r.loc_error), LocResults);
    valid_lr = LocResults(valid_mask);

    if ~isempty(valid_lr)
        all_bands  = [valid_lr.band_id];
        all_errors = [valid_lr.loc_error];
        all_types  = {valid_lr.type_hat};
        unique_bands = unique(all_bands);

        figure('Name', '分频带定位误差', 'Position', [100 50 1400 800]);

        subplot(2, 3, 1);
        grp_data = []; grp_label = [];
        for bi = 1:numel(unique_bands)
            bb = unique_bands(bi);
            errs = all_errors(all_bands == bb);
            grp_data  = [grp_data, errs]; %#ok<AGROW>
            grp_label = [grp_label, bb * ones(1, numel(errs))]; %#ok<AGROW>
        end
        if numel(grp_data) > 1
            boxplot(grp_data, grp_label);
        else
            bar(grp_label, grp_data);
        end
        xlabel('频带'); ylabel('定位误差 (m)');
        title('各频带误差箱线图');
        grid on;

        subplot(2, 3, 2);
        hold on;
        band_colors = [0.2 0.6 1; 0 0.8 0.4; 1 0.7 0; 0.9 0.2 0.2];
        leg_items = {};
        for bi = 1:numel(unique_bands)
            bb = unique_bands(bi);
            errs = sort(all_errors(all_bands == bb));
            n = numel(errs);
            if n == 0, continue; end
            cdf_y = (1:n) / n;
            plot(errs, cdf_y, '-', 'Color', band_colors(bb,:), 'LineWidth', 2);
            leg_items{end+1} = sprintf('Band %d (n=%d, med=%.1fm)', bb, n, median(errs)); %#ok<AGROW>
        end
        hold off;
        xlabel('定位误差 (m)'); ylabel('CDF');
        title('各频带误差累积分布');
        legend(leg_items, 'Location', 'southeast');
        grid on;

        subplot(2, 3, 3);
        rmse_per_band = zeros(1, numel(unique_bands));
        mean_per_band = zeros(1, numel(unique_bands));
        n_per_band    = zeros(1, numel(unique_bands));
        for bi = 1:numel(unique_bands)
            bb = unique_bands(bi);
            errs = all_errors(all_bands == bb);
            rmse_per_band(bi) = sqrt(mean(errs.^2));
            mean_per_band(bi) = mean(errs);
            n_per_band(bi)    = numel(errs);
        end
        bar_data = [mean_per_band; rmse_per_band]';
        bh = bar(unique_bands, bar_data);
        bh(1).FaceColor = [0.3 0.6 0.9];
        bh(2).FaceColor = [0.9 0.3 0.3];
        xlabel('频带'); ylabel('误差 (m)');
        title('各频带 Mean / RMSE');
        legend('Mean', 'RMSE', 'Location', 'best');
        grid on;
        for bi = 1:numel(unique_bands)
            text(unique_bands(bi), rmse_per_band(bi) + 2, ...
                sprintf('n=%d', n_per_band(bi)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end

        subplot(2, 3, 4);
        type_list = {'target', 'opportunistic'};
        type_short = {'target', 'opportunistic'};
        type_colors = [0.9 0.2 0.2; 1 0.7 0];
        hold on;
        legend_h = gobjects(0);
        legend_s = {};
        for ti = 1:numel(type_list)
            tmask = strcmp(all_types, type_list{ti});
            for bi = 1:numel(unique_bands)
                bb = unique_bands(bi);
                bmask = all_bands == bb;
                errs = all_errors(tmask & bmask);
                if isempty(errs), continue; end
                x_pos = bb + (ti - 1.5) * 0.2;
                h = plot(x_pos * ones(size(errs)), errs, 'o', ...
                    'Color', type_colors(ti,:), 'MarkerSize', 5, ...
                    'MarkerFaceColor', type_colors(ti,:));
                plot(x_pos + [-0.1, 0.1], [mean(errs), mean(errs)], '-', ...
                    'Color', type_colors(ti,:)*0.6, 'LineWidth', 2);
                if bi == 1
                    legend_h(end+1) = h; %#ok<AGROW>
                    legend_s{end+1} = type_short{ti}; %#ok<AGROW>
                end
            end
        end
        hold off;
        xlabel('频带'); ylabel('定位误差 (m)');
        title('按源类型分频带误差');
        set(gca, 'XTick', unique_bands);
        if ~isempty(legend_h)
            legend(legend_h, legend_s, 'Location', 'best');
        end
        grid on;

        subplot(2, 3, 5);
        loc_events = EventList(arrayfun(@(ev) ...
            strcmp(ev.route_action, 'localize_only') || ...
            strcmp(ev.route_action, 'localize_then_calibrate'), EventList));
        if numel(loc_events) == numel(valid_lr)
            durations = [loc_events.duration];
            scatter(durations, all_errors, 30, all_bands, 'filled');
            colorbar; clim([0.5 B+0.5]);
            xlabel('事件持续帧数'); ylabel('定位误差 (m)');
            title('误差 vs 事件长度 (颜色=频带)');
            grid on;
        end

        subplot(2, 3, 6);
        if numel(loc_events) == numel(valid_lr)
            powers = [loc_events.power_level_est];
            scatter(powers, all_errors, 30, all_bands, 'filled');
            colorbar; clim([0.5 B+0.5]);
            xlabel('事件功率等级 (dBm)'); ylabel('定位误差 (m)');
            title('误差 vs 功率 (颜色=频带)');
            grid on;
        end

        sgtitle('M4 分频带定位误差统计');

        fprintf('\n--- 分频带定位误差详细统计 ---\n');
        fprintf('  %-6s %-6s %-10s %-10s %-10s %-10s %-10s\n', ...
            'Band', 'N', 'Mean(m)', 'Median(m)', 'RMSE(m)', 'Max(m)', 'Min(m)');
        for bi = 1:numel(unique_bands)
            bb = unique_bands(bi);
            errs = all_errors(all_bands == bb);
            fprintf('  %-6d %-6d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n', ...
                bb, numel(errs), mean(errs), median(errs), ...
                sqrt(mean(errs.^2)), max(errs), min(errs));
        end
        fprintf('  %-6s %-6d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n', ...
            'ALL', numel(all_errors), mean(all_errors), median(all_errors), ...
            sqrt(mean(all_errors.^2)), max(all_errors), min(all_errors));
    end
end

fprintf('\n============================================\n');
fprintf('  仿真完成（新框架：外部标签 + rf_minmax 指纹）\n');
fprintf('============================================\n');


%% ==================== 辅助函数 ====================

function ExternalLabelsRaw = build_labels_from_truth_log(M0Logs, Config)
% BUILD_LABELS_FROM_TRUTH_LOG  从 M0 真值日志提取外部标签（新标签体系）
%
%   使用数字标签 (1/2/3/4)，直接从 source_type 映射。

    LC = source_label_constants();

    ExternalLabelsRaw = struct( ...
        'source_uid',    {}, ...
        'band_id',       {}, ...
        'label',         {}, ...
        'start_frame',   {}, ...
        'end_frame',     {}, ...
        'position_hint', {}, ...
        'template_key',  {}, ...
        'location_prior',{}, ...
        'power_prior',   {}, ...
        'metadata',      {});

    if isempty(M0Logs.TruthLogAll)
        return;
    end

    logs = M0Logs.TruthLogAll;
    keys = {logs.template_key};
    bands = [logs.band_id];
    frames = [logs.frame_id];
    types = {logs.source_type};
    positions = reshape([logs.true_pos_xy], 2, [])';

    inst_ids = [logs.instance_id];
    combo_keys = {};
    for i = 1:numel(keys)
        combo_keys{i} = sprintf('%s__inst%d__b%d', types{i}, inst_ids(i), bands(i)); %#ok<AGROW>
    end
    [unique_combos, ~, ic] = unique(combo_keys);

    for j = 1:numel(unique_combos)
        idx = find(ic == j);
        [seg_frames, ord] = sort(frames(idx));
        idx = idx(ord);
        break_pos = [0, find(diff(seg_frames) > 1), numel(seg_frames)];

        for si = 1:numel(break_pos)-1
            seg_idx = idx((break_pos(si)+1):break_pos(si+1));
            lbl = struct();
            lbl.source_uid = sprintf('%s__seg%d', unique_combos{j}, si);
            lbl.band_id = bands(seg_idx(1));
            lbl.start_frame = min(frames(seg_idx));
            lbl.end_frame = max(frames(seg_idx));
            lbl.position_hint = positions(seg_idx(1), :);
            lbl.template_key = keys{seg_idx(1)};
            lbl.metadata = struct( ...
                'instance_id', inst_ids(seg_idx(1)), ...
                'source_type', types{seg_idx(1)}, ...
                'n_frames', numel(seg_idx));

            src_type = types{seg_idx(1)};
            if LC.type_to_label.isKey(src_type)
                lbl.label = LC.type_to_label(src_type);
            else
                lbl.label = LC.TARGET;
            end

            lbl.location_prior = struct('type', 'none', 'value', []);
            lbl.power_prior    = struct('type', 'none', 'value', []);

            if lbl.label == LC.PERSISTENT_CAL || lbl.label == LC.BROADBAND_CAL
                lbl.location_prior = struct('type', 'exact', 'value', lbl.position_hint);
                lbl.power_prior    = struct('type', 'exact', 'value', logs(seg_idx(1)).tx_power_dBm);
            end

            if isempty(ExternalLabelsRaw)
                ExternalLabelsRaw = lbl;
            else
                ExternalLabelsRaw(end+1) = lbl; %#ok<AGROW>
            end
        end
    end
end


function EventList = build_event_list_from_context(SourceContext, Y_dBm_all, Y_lin_all, Config)
% BUILD_EVENT_LIST_FROM_CONTEXT  从 SourceContext 构建 M4 兼容的 EventList（新标签体系）

    LC = source_label_constants();

    EventList = [];
    if isempty(SourceContext.label_table)
        fprintf('[BuildEvents] 无外部标签，EventList 为空\n');
        return;
    end

    [M_ap, ~, ~] = size(Y_dBm_all);

    for i = 1:numel(SourceContext.label_table)
        lbl = SourceContext.label_table(i);
        ev = struct();
        ev.event_id = i;
        ev.band_id = lbl.band_id;
        ev.t_start = lbl.start_frame;
        ev.t_end = lbl.end_frame;
        ev.time_range = [lbl.start_frame, lbl.end_frame];
        ev.duration = lbl.end_frame - lbl.start_frame + 1;

        b = lbl.band_id;
        ts = lbl.start_frame;
        te = lbl.end_frame;
        ev.obs_segment_dBm = squeeze(Y_dBm_all(:, b, ts:te));
        ev.obs_segment_lin = squeeze(Y_lin_all(:, b, ts:te));
        if ev.duration == 1
            ev.obs_segment_dBm = ev.obs_segment_dBm(:);
            ev.obs_segment_lin = ev.obs_segment_lin(:);
        end

        ev.label = lbl.label;
        if LC.name_map.isKey(lbl.label)
            ev.type_hat = LC.name_map(lbl.label);
        else
            ev.type_hat = 'unknown';
        end
        if LC.route_map.isKey(lbl.label)
            ev.route_action = LC.route_map(lbl.label);
        else
            ev.route_action = 'hold';
        end

        ev.linked_template_key = lbl.template_key;
        ev.source_uid = lbl.source_uid;
        ev.hold_reason = '';
        ev.upgrade_hint = 'none';

        P_bar = mean(ev.obs_segment_dBm, 2);
        ev.power_level_est = mean(P_bar);
        if ev.duration > 1
            ev.power_stability_est = var(mean(ev.obs_segment_dBm, 1));
        else
            ev.power_stability_est = 0;
        end

        ev.score_trusted = 0;
        ev.score_prior_pos = 0;
        ev.score_prior_time = 0;
        ev.score_target = 0;
        ev.band_coverage_vec = zeros(1, Config.m0.num_bands);
        ev.band_coverage_vec(b) = 1;
        ev.n_valid_ap = M_ap;

        if isempty(EventList)
            EventList = ev;
        else
            EventList(end+1) = ev; %#ok<AGROW>
        end
    end

    fprintf('[BuildEvents] 从 %d 条标签构建 %d 个事件\n', ...
        numel(SourceContext.label_table), numel(EventList));
end


function h = plot_unique_sources(all_xy, all_ids, mask, marker, color, sz)
    h = [];
    if ~any(mask), return; end
    ids = unique(all_ids(mask));
    for ii = 1:numel(ids)
        idx_first = find(all_ids == ids(ii) & mask, 1);
        xy = all_xy(idx_first, :);
        if ii == 1
            h = plot(xy(1), xy(2), marker, 'Color', color * 0.8, ...
                'MarkerSize', sz, 'MarkerFaceColor', color, 'LineWidth', 1.2);
        else
            plot(xy(1), xy(2), marker, 'Color', color * 0.8, ...
                'MarkerSize', sz, 'MarkerFaceColor', color, 'LineWidth', 1.2);
        end
    end
end
