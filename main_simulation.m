%% MAIN_SIMULATION  主仿真入口 — M0 源调度 + M1 信道观测 + M2.5 指纹库
%
%  功能：
%    1. M0 初始化 + M1 初始化 + M2.5 指纹库构建
%    2. 主循环：逐帧调度源 → 生成 RSS 观测
%    3. 硬约束验证
%    4. 统计输出
%    5. 可视化（场景布局、频带占用时间线、RSS 热力图、RSS vs 距离、指纹热力图）
%
%  后续模块接入点：在主循环内 step_m1 之后添加 M3/M4 等

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  频谱指纹匹配机遇式增强定位 — 主仿真\n');
fprintf('============================================\n\n');

%% ========== 第一阶段：初始化 ==========

% --- M0 初始化 ---
% --- 场景配置（可调） ---
sim_override.ap.num_x     = 8;   % 8 x 8 = 64 个 AP
sim_override.ap.num_y     = 4;
sim_override.fp.grid_step = 1;   % 1m 网格

[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);

% --- M1 初始化 ---
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

% --- M2.5 指纹库与模板库 ---
[SpatialFP, SignatureLib] = init_m25_single_source_fp( ...
    APs, Bands, GridValid, Config, SourceTemplates);

% --- 全局参数 ---
T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;
dt = Config.m0.dt;

fprintf('\n--- 仿真参数 ---\n');
fprintf('  区域: [%.0f, %.0f] x [%.0f, %.0f] m\n', Config.area.x_range, Config.area.y_range);
fprintf('  AP: %d 个 | 有效网格: %d 点\n', M, GridValid.Nvalid);
fprintf('  频带: %d 个 | 总帧数: %d | dt=%.1f s\n', B, T, dt);
for b = 1:B
    fprintf('  Band %d: %s, fc=%.1f MHz, model=%s\n', ...
        b, Bands.name{b}, Bands.fc_Hz(b)/1e6, Bands.model{b});
end

%% ========== 第二阶段：主循环 ==========

fprintf('\n--- 主循环开始 (%d 帧) ---\n', T);

% 预分配存储
FrameStates    = cell(1, T);
ObsFrames      = cell(1, T);
Y_dBm_all      = zeros(M, B, T);
Y_lin_all       = zeros(M, B, T);
occupancy_matrix = zeros(T, B);   % 0=空, 1=trusted, 2=prior_pos, 3=prior_time, 4=target

tic;
for t = 1:T
    % --- M0：源调度 ---
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);

    % --- M1：信道观测 ---
    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);

    % --- 后续模块（M3/M4 在主循环后以事件级处理） ---
    % [TODO] M5/M6: 误差统计与标校

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
                case 'trusted_fixed',   occupancy_matrix(t, b) = 1;
                case 'prior'
                    if strcmp(apb.source_subtype, 'prior_pos_known')
                        occupancy_matrix(t, b) = 2;
                    else
                        occupancy_matrix(t, b) = 3;
                    end
                case 'ordinary_target', occupancy_matrix(t, b) = 4;
            end
        end
    end
end
elapsed = toc;
fprintf('主循环完成: %d 帧, 耗时 %.2f s (%.1f 帧/s)\n', T, elapsed, T/elapsed);

%% ========== 第 2.5 阶段：M3 事件分类 + M4 WKNN 定位 ==========

% --- M3：单源事件分类 ---
[EventList, GroupListM3] = run_m3_event_classifier_single_source( ...
    Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config);

% --- M4：WKNN 定位（只对需定位的事件） ---
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
type_names = {'空闲', 'trusted_fixed', 'prior_pos_known', 'prior_time_known', 'ordinary_target'};
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
plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

if ~isempty(M0Logs.TruthLogAll)
    all_types = {M0Logs.TruthLogAll.source_type};
    all_xy    = reshape([M0Logs.TruthLogAll.true_pos_xy], 2, [])';
    all_ids   = [M0Logs.TruthLogAll.instance_id];
    all_keys  = {M0Logs.TruthLogAll.template_key};

    % trusted_fixed (蓝色五角星)
    mask_tr = strcmp(all_types, 'trusted_fixed');
    h_tr = plot_unique_sources(all_xy, all_ids, mask_tr, 'p', [0.2 0.5 1], 16);

    % prior_pos_known (绿色方块)
    mask_pp = strncmp(all_keys, 'prior_pos_', 10);
    h_pp = plot_unique_sources(all_xy, all_ids, mask_pp, 's', [0 0.8 0.4], 10);

    % prior_time_known (橙色菱形)
    mask_pt = strncmp(all_keys, 'prior_time_', 11);
    h_pt = plot_unique_sources(all_xy, all_ids, mask_pt, 'd', [1 0.7 0], 9);

    % ordinary_target (深色圆点)
    mask_tg = strcmp(all_types, 'ordinary_target');
    h_tg = plot_unique_sources(all_xy, all_ids, mask_tg, 'o', [0.1 0.2 0.2], 7);
end

xlabel('X (m)'); ylabel('Y (m)');
title('场景布局 — 各类源累积出现位置');
axis equal; grid on;
xlim(Config.area.x_range); ylim(Config.area.y_range);

% 图例
leg_h = {}; leg_s = {};
leg_h{end+1} = plot(NaN, NaN, '.', 'Color', [0.85 0.85 0.85], 'MarkerSize', 10);
leg_s{end+1} = '有效网格点';
leg_h{end+1} = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
leg_s{end+1} = 'AP';
if ~isempty(h_tr), leg_h{end+1} = h_tr; leg_s{end+1} = 'trusted\_fixed'; end
if ~isempty(h_pp), leg_h{end+1} = h_pp; leg_s{end+1} = 'prior\_pos\_known'; end
if ~isempty(h_pt), leg_h{end+1} = h_pt; leg_s{end+1} = 'prior\_time\_known'; end
if ~isempty(h_tg), leg_h{end+1} = h_tg; leg_s{end+1} = 'ordinary\_target'; end
legend([leg_h{:}], leg_s, 'Location', 'best');

% ===== 图4：最后有源帧 RSS vs 距离 =====
figure('Name', 'RSS vs 距离', 'Position', [700 700 1000 400]);
plot_count = 0;
for b = 1:B
    src_frames = find(occupancy_matrix(:, b) > 0);
    if isempty(src_frames), continue; end
    t_last = src_frames(end);

    src_xy = FrameStates{t_last}.active_per_band(b).true_pos_xy;
    dist_m = sqrt(sum((APs.pos_xy - src_xy).^2, 2));

    plot_count = plot_count + 1;
    subplot(1, B, b);
    scatter(dist_m, Y_dBm_all(:, b, t_last), 40, 'filled');
    xlabel('距离 (m)'); ylabel('RSS (dBm)');
    title(sprintf('Band %d 帧%d', b, t_last));
    grid on;
end
sgtitle('RSS vs 源-AP 距离');

% ===== 图5：SpatialFP 指纹热力图 =====
% 含义：在每个网格点放置参考信源(0dBm)，16个AP接收的平均RSS
% mean_dBm(g) 反映该位置信源被 AP 阵列整体感知的强弱
figure('Name', '指纹库热力图', 'Position', [100 100 1200 700]);
for b = 1:B
    subplot(2, 2, b);
    scatter(SpatialFP.grid_xy(:,1), SpatialFP.grid_xy(:,2), ...
        12, SpatialFP.band(b).mean_dBm, 'filled');
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    colorbar;
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('Band %d (%s) 各点平均 RSS [dBm]', b, Bands.name{b}));
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);
end
sgtitle(sprintf('M2.5 指纹库: 网格点放信源(ref=%d dBm), 16AP 平均接收 RSS', ...
    SpatialFP.ref_power_dBm(1)));

% ===== 图6：M3 事件分类时间线（与图1对齐的 imagesc 风格） =====
if ~isempty(EventList)
    figure('Name', 'M3 事件分类', 'Position', [50 50 1200 350]);

    % 构建 T x B 矩阵：0=空闲, 1=trusted, 2=prior_pos, 3=prior_time, 4=target
    m3_matrix = zeros(T, B);
    for e = 1:numel(EventList)
        ev = EventList(e);
        frames = ev.t_start:ev.t_end;
        switch ev.type_hat
            case 'trusted_fixed',    val = 1;
            case 'prior_pos_known',  val = 2;
            case 'prior_time_known', val = 3;
            case 'ordinary_target',  val = 4;
            otherwise,               val = 0;
        end
        m3_matrix(frames, ev.band_id) = val;
    end

    imagesc(1:T, 1:B, m3_matrix');
    colormap(gca, cmap);
    clim([0 4]);
    cbar6 = colorbar;
    cbar6.Ticks = [0, 1, 2, 3, 4];
    cbar6.TickLabels = type_names;
    xlabel('帧号'); ylabel('频带');
    title('M3 事件分类时间线');
    set(gca, 'YTick', 1:B);
    grid on;
end

% ===== 图7：M4 定位结果散点图 =====
if ~isempty(LocResults)
    figure('Name', 'M4 定位结果', 'Position', [100 100 800 700]);

    % 背景：有效网格 + AP
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    % 画每个定位结果
    for k = 1:numel(LocResults)
        lr = LocResults(k);
        if isempty(lr.true_pos_xy) || all(lr.true_pos_xy == 0)
            continue;
        end

        switch lr.type_hat
            case 'ordinary_target',  clr = [0.9 0.2 0.2];
            case 'prior_time_known', clr = [1 0.7 0];
            otherwise,               clr = [0.5 0.5 0.5];
        end

        % 真值位置
        plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', ...
            'Color', clr, 'MarkerSize', 8, 'MarkerFaceColor', clr);
        % 估计位置
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
            'Color', clr*0.6, 'MarkerSize', 10, 'LineWidth', 2);
        % 连线
        plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
             [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
             '-', 'Color', [clr, 0.4], 'LineWidth', 1);
    end

    xlabel('X (m)'); ylabel('Y (m)');
    title('M4 WKNN 定位结果 (o=真值, x=估计)');
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);

    % 图例（在 hold on 状态下画不可见标记）
    leg_h = gobjects(5, 1);
    leg_h(1) = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    leg_h(2) = plot(NaN, NaN, 'o', 'Color', [0.9 0.2 0.2], 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.2 0.2]);
    leg_h(3) = plot(NaN, NaN, 'x', 'Color', [0.9 0.2 0.2]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    leg_h(4) = plot(NaN, NaN, 'o', 'Color', [1 0.7 0], 'MarkerSize', 8, 'MarkerFaceColor', [1 0.7 0]);
    leg_h(5) = plot(NaN, NaN, 'x', 'Color', [1 0.7 0]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    legend(leg_h, {'AP', 'target 真值', 'target 估计', 'prior\_time 真值', 'prior\_time 估计'}, ...
        'Location', 'best');
end

% ===== 图8：分频带定位误差统计 =====
if ~isempty(LocResults)
    % 收集有真值的定位结果
    valid_mask = ~arrayfun(@(r) isnan(r.loc_error), LocResults);
    valid_lr = LocResults(valid_mask);

    if ~isempty(valid_lr)
        all_bands  = [valid_lr.band_id];
        all_errors = [valid_lr.loc_error];
        all_types  = {valid_lr.type_hat};
        unique_bands = unique(all_bands);

        figure('Name', '分频带定位误差', 'Position', [100 50 1400 800]);

        % --- 子图1：各频带误差箱线图 ---
        subplot(2, 3, 1);
        grp_data = [];
        grp_label = [];
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

        % --- 子图2：各频带误差 CDF ---
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

        % --- 子图3：各频带 RMSE 柱状图 ---
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
        % 标注样本数
        for bi = 1:numel(unique_bands)
            text(unique_bands(bi), rmse_per_band(bi) + 2, ...
                sprintf('n=%d', n_per_band(bi)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end

        % --- 子图4：按源类型分频带误差 ---
        subplot(2, 3, 4);
        type_list = {'ordinary_target', 'prior_time_known'};
        type_short = {'target', 'prior\_time'};
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
                % 画均值线
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

        % --- 子图5：误差 vs 事件持续帧数 ---
        subplot(2, 3, 5);
        valid_events = EventList(arrayfun(@(ev) ...
            strcmp(ev.route_action, 'localize_only') || ...
            strcmp(ev.route_action, 'localize_then_calibrate'), EventList));
        if numel(valid_events) == numel(valid_lr)
            durations = [valid_events.duration];
            scatter(durations, all_errors, 30, all_bands, 'filled');
            colorbar; clim([0.5 B+0.5]);
            xlabel('事件持续帧数'); ylabel('定位误差 (m)');
            title('误差 vs 事件长度 (颜色=频带)');
            grid on;
        end

        % --- 子图6：误差 vs 估计功率等级 ---
        subplot(2, 3, 6);
        if numel(valid_events) == numel(valid_lr)
            powers = [valid_events.power_level_est];
            scatter(powers, all_errors, 30, all_bands, 'filled');
            colorbar; clim([0.5 B+0.5]);
            xlabel('事件功率等级 (dBm)'); ylabel('定位误差 (m)');
            title('误差 vs 功率 (颜色=频带)');
            grid on;
        end

        sgtitle('M4 分频带定位误差统计');

        % --- 终端打印详细统计 ---
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
fprintf('  仿真完成\n');
fprintf('============================================\n');


%% ==================== 辅助函数 ====================

function h = plot_unique_sources(all_xy, all_ids, mask, marker, color, sz)
% PLOT_UNIQUE_SOURCES  按实例去重后绘制源位置
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
