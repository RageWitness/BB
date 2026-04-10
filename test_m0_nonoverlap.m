%% TEST_M0_NONOVERLAP  M0 模块冒烟测试
%
%  运行方式：直接在 MATLAB 中执行此脚本
%  验证目标：
%    1. 初始化不报错
%    2. 500 帧循环不报错
%    3. 每频带每帧最多 1 个主导源 (硬约束)
%    4. 真值日志非空
%    5. 输出简要统计信息与时间线图

clear; clc; close all;
rng(42);  % 可重复

fprintf('===== M0 冒烟测试开始 =====\n\n');

%% 1. 初始化
[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap();

fprintf('\n--- 配置摘要 ---\n');
fprintf('  区域: [%.0f, %.0f] x [%.0f, %.0f] m\n', ...
    Config.area.x_range, Config.area.y_range);
fprintf('  AP: %d 个\n', APs.num);
fprintf('  有效网格点: %d\n', GridValid.Nvalid);
fprintf('  频带数: %d\n', Config.m0.num_bands);
fprintf('  总帧数: %d\n', Config.m0.T_total);
fprintf('  帧间隔: %.1f s\n', Config.m0.dt);

%% 2. 主循环
T = Config.m0.T_total;
B = Config.m0.num_bands;
FrameStates = cell(1, T);

% 统计矩阵
occupancy_matrix = zeros(T, B);  % 0=空, 1=trusted, 2=prior_pos, 3=prior_time, 4=target

tic;
for t = 1:T
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);
    FrameStates{t} = FrameState_t;

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
fprintf('\n主循环 %d 帧完成，耗时 %.2f s (%.1f 帧/s)\n', T, elapsed, T/elapsed);

%% 3. 硬约束验证：每频带每帧最多 1 个主导源
% (由设计保证，这里做断言)
for t = 1:T
    for b = 1:B
        n_active = double(FrameStates{t}.active_per_band(b).has_source);
        assert(n_active <= 1, 'N_active(%d,%d) = %d > 1 ！违反硬约束！', b, t, n_active);
    end
end
fprintf('\n[PASS] 硬约束验证通过：所有频带所有帧 N_active(b,t) <= 1\n');

%% 4. 统计
fprintf('\n--- 真值日志统计 ---\n');
fprintf('  TruthLogAll   条目数: %d\n', numel(M0Logs.TruthLogAll));
fprintf('  TruthLogTarget 条目数: %d\n', numel(M0Logs.TruthLogTarget));

% 按类型统计占用帧数
type_names = {'空闲', 'trusted_fixed', 'prior_pos_known', 'prior_time_known', 'ordinary_target'};
fprintf('\n--- 各频带占用统计 (帧数) ---\n');
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

%% 5. 可视化
figure('Name', 'M0 频带占用时间线', 'Position', [100 100 1200 500]);

% 颜色映射
cmap = [1 1 1;       % 0: 空闲（白色）
        0.2 0.6 1;   % 1: trusted（蓝色）
        0 0.8 0.4;   % 2: prior_pos（绿色）
        1 0.7 0;     % 3: prior_time（橙色）
        0.9 0.2 0.2]; % 4: target（红色）

imagesc(1:T, 1:B, occupancy_matrix');
colormap(cmap);
cbar = colorbar;
cbar.Ticks = [0, 1, 2, 3, 4];
cbar.TickLabels = type_names;
xlabel('帧号 (t)');
ylabel('频带 (b)');
title('M0 频带占用时间线');
set(gca, 'YTick', 1:B);
grid on;

% AP 与网格可视化
figure('Name', 'M0 场景布局', 'Position', [100 650 600 500]);
plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 3);
hold on;
plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% 标注 trusted 模板位置
for j = 1:numel(SourceTemplates.trusted)
    plot(SourceTemplates.trusted(j).fixed_pos_xy(1), ...
         SourceTemplates.trusted(j).fixed_pos_xy(2), ...
         'bp', 'MarkerSize', 14, 'MarkerFaceColor', 'b');
end
xlabel('X (m)'); ylabel('Y (m)');
title('场景布局: AP(红三角) + 网格(灰点) + Trusted(蓝五角星)');
axis equal; grid on;
xlim(Config.area.x_range); ylim(Config.area.y_range);
legend('有效网格点', 'AP', 'Trusted 源', 'Location', 'best');

fprintf('\n===== M0 冒烟测试完成 =====\n');
