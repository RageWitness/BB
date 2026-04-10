%% TEST_M1_CHANNEL  M1 信道模块冒烟测试
%
%  运行方式：直接在 MATLAB 中执行此脚本
%  前置条件：需要 Communications Toolbox（awgn 函数）
%  验证目标：
%    1. M0 + M1 联合初始化不报错
%    2. 100 帧循环不报错
%    3. Y_dBm 维度正确 = 16 x 4
%    4. 有源频带信号显著高于噪声底
%    5. 可视化：RSS 热力图 + RSS vs 距离散点图

clear; clc; close all;
rng(42);

fprintf('===== M1 信道模块冒烟测试 =====\n\n');

%% 1. 初始化 M0
[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap();

%% 2. 初始化 M1
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

%% 3. 主循环
T_test = 100;
M = APs.num;
B = Bands.B;

Y_dBm_all = zeros(M, B, T_test);
has_source_all = false(B, T_test);

fprintf('\n--- 开始 %d 帧循环 ---\n', T_test);
tic;
for t = 1:T_test
    % M0 步进
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);

    % M1 步进
    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);

    % 收集结果
    Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;
    for b = 1:B
        has_source_all(b, t) = FrameState_t.active_per_band(b).has_source;
    end
end
elapsed = toc;
fprintf('主循环完成，耗时 %.2f s (%.1f 帧/s)\n', elapsed, T_test/elapsed);

%% 4. 基本验证
fprintf('\n--- 基本验证 ---\n');
assert(all(size(ObsFrame.Y_dBm) == [M, B]), 'Y_dBm 维度错误！');
fprintf('[PASS] Y_dBm 维度: %d x %d\n', size(ObsFrame.Y_dBm));

% 检查有源帧的信号是否高于噪声底
for b = 1:B
    src_frames = find(has_source_all(b, :));
    if ~isempty(src_frames)
        mean_src = mean(mean(Y_dBm_all(:, b, src_frames)));
        nosrc_frames = find(~has_source_all(b, :));
        if ~isempty(nosrc_frames)
            mean_nosrc = mean(mean(Y_dBm_all(:, b, nosrc_frames)));
            fprintf('  Band %d (%s): 有源均值=%.1f dBm, 无源均值=%.1f dBm, 差=%.1f dB\n', ...
                b, Bands.name{b}, mean_src, mean_nosrc, mean_src - mean_nosrc);
        else
            fprintf('  Band %d (%s): 有源均值=%.1f dBm (全部帧有源)\n', ...
                b, Bands.name{b}, mean_src);
        end
    else
        fprintf('  Band %d (%s): 无有源帧\n', b, Bands.name{b});
    end
end

%% 5. 统计
fprintf('\n--- 源出现统计 ---\n');
for b = 1:B
    n_src = sum(has_source_all(b, :));
    fprintf('  Band %d: %d/%d 帧有源 (%.0f%%)\n', b, n_src, T_test, 100*n_src/T_test);
end

%% 6. 可视化







% --- 图1：RSS 热力图 (AP x 帧) ---
figure('Name', 'M1 RSS 热力图', 'Position', [50 50 1200 700]);
for b = 1:B
    subplot(2, 2, b);
    imagesc(1:T_test, 1:M, squeeze(Y_dBm_all(:, b, :)));
    colorbar;
    xlabel('帧号'); ylabel('AP');
    title(sprintf('Band %d (%s) RSS [dBm]', b, Bands.name{b}));
    set(gca, 'YTick', 1:4:M);
end
sgtitle('M1 RSS 热力图 (AP x 帧)');

% --- 图2：最后一个有源帧 RSS vs 距离 ---
figure('Name', 'M1 RSS vs 距离', 'Position', [50 800 1200 500]);
for b = 1:B
    src_frames = find(has_source_all(b, :));
    if isempty(src_frames), continue; end
    t_last = src_frames(end);

    % 重新跑该帧获取 meta
    [~, FS_tmp, ~] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t_last, M0Logs);
    src_xy = FS_tmp.active_per_band(b).true_pos_xy;
    dist_m = sqrt(sum((APs.pos_xy - src_xy).^2, 2));

    subplot(1, B, b);
    scatter(dist_m, Y_dBm_all(:, b, t_last), 40, 'filled');
    hold on;
    xlabel('距离 (m)'); ylabel('RSS (dBm)');
    title(sprintf('Band %d (%s) 帧%d', b, Bands.name{b}, t_last));
    grid on;
end
sgtitle('RSS vs 源-AP 距离');

fprintf('\n===== M1 冒烟测试完成 =====\n');



