%% MAIN_SIMULATION  主仿真入口 — M0 源调度 + M1 信道观测
%
%  功能：
%    1. M0 初始化 + M1 初始化
%    2. 主循环：逐帧调度源 → 生成 RSS 观测
%    3. 硬约束验证
%    4. 统计输出
%    5. 可视化（场景布局、频带占用时间线、RSS 热力图、RSS vs 距离）
%
%  后续模块接入点：在主循环内 step_m1 之后添加 M2.5/M3/M4 等

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  频谱指纹匹配机遇式增强定位 — 主仿真\n');
fprintf('============================================\n\n');

%% ========== 第一阶段：初始化 ==========

% --- M0 初始化 ---
[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap();

% --- M1 初始化 ---
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

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
occupancy_matrix = zeros(T, B);   % 0=空, 1=trusted, 2=prior_pos, 3=prior_time, 4=target

tic;
for t = 1:T
    % --- M0：源调度 ---
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);

    % --- M1：信道观测 ---
    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);

    % --- 后续模块接入点 ---
    % [TODO] M2.5: 指纹库查询
    % [TODO] M3:   单源事件检测
    % [TODO] M4:   WKNN 定位
    % [TODO] M5/M6: 误差统计与标校

    % --- 存储 ---
    FrameStates{t} = FrameState_t;
    ObsFrames{t}   = ObsFrame;
    Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;

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
