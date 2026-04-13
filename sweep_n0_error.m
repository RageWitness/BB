%% SWEEP_N0_ERROR  扫描 AWGN 噪声功率谱密度 n0，观察定位误差变化
%
%  独立脚本，无需先运行 main_simulation.m
%  流程：
%    1. 一次性初始化（M0 / M1 / M2.5）
%    2. 运行 M0 循环，收集全部 FrameStates（与 n0 无关）
%    3. 对每个 n0 值：恢复 RNG → 重跑 M1 观测 → 跳过 M3 直接 WKNN 定位
%    4. 画 n0-误差曲线 + 分频带对比图

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  n0 噪声功率扫描 — 定位误差分析\n');
fprintf('============================================\n\n');

%% ========== 1. 一次性初始化 ==========

% --- 场景配置（可调） ---
ap_override.ap.num_x   = 8;   
ap_override.ap.num_y   = 4;
ap_override.fp.grid_step = 1;  % 1m 网格（~90000 有效点）

[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(ap_override);
[Bands, ~, Config] = init_m1_channel(Config, APs);
[SpatialFP, ~] = init_m25_single_source_fp(APs, Bands, GridValid, Config, SourceTemplates);

T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;
K = 3;  % WKNN 邻居数

%% ========== 2. 运行 M0 循环（与 n0 无关） ==========

fprintf('\n--- 运行 M0 源调度 (%d 帧) ---\n', T);
FrameStates = cell(1, T);
for t = 1:T
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);
    FrameStates{t} = FrameState_t;
end
fprintf('M0 完成\n');

% 保存 RNG 状态 — M1 循环起点
rng_state_after_m0 = rng;

%% ========== 3. n0 扫描 ==========

n0_sweep = -180 : 4 : -50;    % dBm/Hz: 覆盖全频带 SNR 从 >100dB 到负值
N_sweep  = numel(n0_sweep);

% 预分配：总体统计
mean_err   = zeros(N_sweep, 1);
rmse_err   = zeros(N_sweep, 1);
median_err = zeros(N_sweep, 1);
p90_err    = zeros(N_sweep, 1);
n_samp     = zeros(N_sweep, 1);

% 预分配：分频带统计
mean_err_band = zeros(N_sweep, B);
rmse_err_band = zeros(N_sweep, B);
n_samp_band   = zeros(N_sweep, B);

fprintf('\n--- 开始 n0 扫描 (%d 个点) ---\n', N_sweep);
sweep_tic = tic;

for si = 1:N_sweep
    n0_val = n0_sweep(si);

    % (a) 恢复 RNG，使阴影衰落实现一致
    rng(rng_state_after_m0);

    % (b) 重建干净的 ChannelState，覆写 n0
    ChannelState = init_m1_channel_state(Config, APs, Bands);
    ChannelState.noise_n0_dBmHz = n0_val * ones(1, B);

    % (c) M1 观测循环
    Y_lin_all = zeros(M, B, T);
    for t = 1:T
        [ChannelState, ObsFrame] = step_m1_generate_obs( ...
            FrameStates{t}, APs, Bands, ChannelState, Config, t);
        Y_lin_all(:, :, t) = ObsFrame.Y_lin;
    end

    % (d) 跳过 M3，直接 WKNN 定位
    errors_all = [];
    bands_all  = [];
    for t = 1:T
        for b = 1:B
            apb = FrameStates{t}.active_per_band(b);
            if ~apb.has_source, continue; end
            if ~strcmp(apb.source_type, 'ordinary_target'), continue; end

            F_obs = Y_lin_all(:, b, t);
            F_shape = F_obs / (sum(F_obs) + 1e-30);
            [D_vec, ~, ~, ~] = compute_wknn_distance_shape_scale( ...
                F_obs, F_shape, SpatialFP.band(b));
            [ni, ~, nw] = select_knn_neighbors_m4(D_vec, K);
            ep = estimate_position_wknn_m4(SpatialFP.grid_xy, ni, nw);

            err = norm(ep - apb.true_pos_xy);
            errors_all(end+1) = err; 
            bands_all(end+1)  = b;   
        end
    end

    % (e) 记录统计
    if ~isempty(errors_all)
        mean_err(si)   = mean(errors_all);
        rmse_err(si)   = sqrt(mean(errors_all.^2));
        median_err(si) = median(errors_all);
        p90_err(si)    = prctile(errors_all, 90);
        n_samp(si)     = numel(errors_all);

        for b = 1:B
            mask = (bands_all == b);
            if any(mask)
                eb = errors_all(mask);
                mean_err_band(si, b) = mean(eb);
                rmse_err_band(si, b) = sqrt(mean(eb.^2));
                n_samp_band(si, b)   = sum(mask);
            end
        end
    end

    fprintf('  [%2d/%d] n0 = %4.0f dBm/Hz | Mean=%.1f m | RMSE=%.1f m | n=%d\n', ...
        si, N_sweep, n0_val, mean_err(si), rmse_err(si), n_samp(si));
end

fprintf('扫描完成, 耗时 %.1f s\n', toc(sweep_tic));

%% ========== 4. 可视化 ==========

% 预计算各频带的近似 SNR 矩阵 (N_sweep x B)
snr_approx_all = zeros(N_sweep, B);
for b = 1:B
    N_dBm_sweep = n0_sweep + 10*log10(Bands.bw_Hz(b));  % 噪声总功率
    mean_rx_dBm = mean(SpatialFP.band(b).mean_dBm);     % 库平均接收功率
    snr_approx_all(:, b) = mean_rx_dBm - N_dBm_sweep';
end

% 预计算各频带的带宽标签字符串
bw_labels = cell(1, B);
for b = 1:B
    hz = Bands.bw_Hz(b);
    if hz >= 1e9
        bw_labels{b} = sprintf('%.1f GHz', hz / 1e9);
    elseif hz >= 1e6
        bw_labels{b} = sprintf('%.1f MHz', hz / 1e6);
    elseif hz >= 1e3
        bw_labels{b} = sprintf('%.1f kHz', hz / 1e3);
    else
        bw_labels{b} = sprintf('%.0f Hz', hz);
    end
end

band_colors = lines(B);

% --- 图1：n0 vs 总体误差 ---
figure('Name', 'n0 vs Localization Error', 'Position', [50 50 900 550]);
plot(n0_sweep, mean_err,   'o-',  'LineWidth', 2, 'DisplayName', 'Mean');
hold on;
plot(n0_sweep, rmse_err,   's-',  'LineWidth', 2, 'DisplayName', 'RMSE');
plot(n0_sweep, median_err, 'd-',  'LineWidth', 2, 'DisplayName', 'Median');
plot(n0_sweep, p90_err,    '^-',  'LineWidth', 2, 'DisplayName', '90th percentile');
xline(-174, 'k--', 'n_0 = -174 (thermal)', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom');
hold off;
xlabel('n_0 (dBm/Hz)'); ylabel('定位误差 (m)');
title('AWGN 噪声功率谱密度 vs 直接 WKNN 定位误差');
legend('Location', 'northwest'); grid on;

% --- 图2：n0 vs 分频带 Mean Error ---
figure('Name', 'Per-band Mean Error vs n0', 'Position', [100 100 900 550]);
hold on;
for b = 1:B
    valid = n_samp_band(:, b) > 0;
    plot(n0_sweep(valid), mean_err_band(valid, b), 'o-', ...
        'Color', band_colors(b,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Band %d (%s)', b, Bands.name{b}));
end
xline(-174, 'k--', 'thermal', 'LineWidth', 1.5);
hold off;
xlabel('n_0 (dBm/Hz)'); ylabel('Mean Error (m)');
title('分频带 Mean Error vs n_0');
legend('Location', 'northwest'); grid on;

% --- 图3：SNR vs 分频带 RMSE（以 SNR 为横轴） ---
figure('Name', 'SNR vs RMSE per band', 'Position', [150 150 900 550]);
hold on;
for b = 1:B
    valid = n_samp_band(:, b) > 0;
    plot(snr_approx_all(valid, b), rmse_err_band(valid, b), 's-', ...
        'Color', band_colors(b,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Band %d (%s, BW=%s)', ...
            b, Bands.name{b}, bw_labels{b}));
end
xline(0, 'r--', 'SNR=0 dB', 'LineWidth', 1.5);
xline(10, 'k:', 'SNR=10 dB', 'LineWidth', 1);
hold off;
xlabel('近似 SNR (dB)'); ylabel('RMSE (m)');
title('各频带 RMSE vs 近似 SNR');
legend('Location', 'northeast'); grid on;
set(gca, 'XDir', 'reverse');  % SNR 高在左

% --- 图4：SNR vs 分频带 Mean Error（以 SNR 为横轴） ---
figure('Name', 'SNR vs Mean Error per band', 'Position', [200 200 900 550]);
hold on;
for b = 1:B
    valid = n_samp_band(:, b) > 0;
    plot(snr_approx_all(valid, b), mean_err_band(valid, b), 'o-', ...
        'Color', band_colors(b,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Band %d (%s)', b, Bands.name{b}));
end
xline(0, 'r--', 'SNR=0 dB', 'LineWidth', 1.5);
xline(10, 'k:', 'SNR=10 dB', 'LineWidth', 1);
hold off;
xlabel('近似 SNR (dB)'); ylabel('Mean Error (m)');
title('各频带 Mean Error vs 近似 SNR');
legend('Location', 'northeast'); grid on;
set(gca, 'XDir', 'reverse');

% --- 图5：分频带双轴图 — 误差 + SNR vs n0 ---
figure('Name', 'Error & SNR vs n0', 'Position', [250 100 1400 500]);
for b = 1:B
    subplot(1, B, b);
    valid = n_samp_band(:, b) > 0;

    % 左轴：误差
    yyaxis left;
    plot(n0_sweep(valid), mean_err_band(valid, b), 'o-', 'LineWidth', 1.5);
    ylabel('Mean Error (m)');

    % 右轴：SNR
    yyaxis right;
    plot(n0_sweep, snr_approx_all(:, b), 's--', 'LineWidth', 1.5);
    ylabel('Approx SNR (dB)');
    yline(0, 'r:', 'SNR=0');

    xlabel('n_0 (dBm/Hz)');
    title(sprintf('Band %d (%s, BW=%s)', b, Bands.name{b}, bw_labels{b}));
    grid on;
end
sgtitle('误差 & 近似 SNR vs n_0');

%% ========== 5. 汇总表 ==========

fprintf('\n--- 扫描结果汇总 ---\n');
fprintf('%-14s %-10s %-10s %-10s %-10s %-6s\n', ...
    'n0(dBm/Hz)', 'Mean(m)', 'RMSE(m)', 'Median(m)', 'P90(m)', 'N');
for si = 1:N_sweep
    fprintf('%-14.0f %-10.2f %-10.2f %-10.2f %-10.2f %-6d\n', ...
        n0_sweep(si), mean_err(si), rmse_err(si), ...
        median_err(si), p90_err(si), n_samp(si));
end
fprintf('\n');

