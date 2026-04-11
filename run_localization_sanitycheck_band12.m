%% RUN_LOCALIZATION_SANITYCHECK_BAND12  Band 1/2 分阶段定位验证
%
%  三组实验，逐步增加难度：
%    A: sigma_dB=0, 低噪声, 固定 0 dBm → 验证匹配逻辑本身
%    B: 恢复 lognormal 阴影, 固定 0 dBm → 验证库/观测是否匹配
%    C: 恢复阴影 + 随机发射功率 + 约束 alpha → 验证未知功率匹配
%
%  只用 Band 1/2（lognormal 频带），跳过 M3 直接单帧 WKNN。

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  Band 1/2 分阶段定位 Sanity Check\n');
fprintf('============================================\n\n');

%% ========== 公共初始化 ==========

base_override.ap.num_x     = 4;
base_override.ap.num_y     = 4;
base_override.fp.grid_step = 1;

% 只用 Band 1/2
test_bands = [1, 2];
K = 5;

% 三组实验配置
experiments = struct();

% --- 实验 A：无阴影、低噪声、固定 0 dBm ---
experiments(1).name     = 'A: 无阴影 + 低噪声 + 0dBm';
experiments(1).sigma_dB = [0, 0, 8, 8];       % Band 1/2 阴影关闭
experiments(1).n0_dBmHz = -200;                % 极低噪声
experiments(1).power_range = [0, 0; 0, 0; -12, -3; -12, -3];  % Band 1/2 固定 0 dBm

% --- 实验 B：恢复阴影、低噪声、固定 0 dBm ---
experiments(2).name     = 'B: 有阴影 + 低噪声 + 0dBm';
experiments(2).sigma_dB = [6, 6, 8, 8];        % 恢复默认阴影
experiments(2).n0_dBmHz = -200;
experiments(2).power_range = [0, 0; 0, 0; -12, -3; -12, -3];

% --- 实验 C：阴影 + 正常噪声 + 随机功率 + 约束 alpha ---
experiments(3).name     = 'C: 有阴影 + 正常噪声 + 随机功率';
experiments(3).sigma_dB = [6, 6, 8, 8];
experiments(3).n0_dBmHz = -174;                % 正常热噪声
experiments(3).power_range = [-15, -5; -15, -5; -12, -3; -12, -3];

% 存储结果
all_results = cell(1, 3);

for ei = 1:3
    exp = experiments(ei);
    fprintf('\n====== 实验 %s ======\n', exp.name);

    % ---- 构建 override ----
    ov = base_override;
    ov.m1.channel.lognormal.sigma_dB = exp.sigma_dB;
    ov.m0.target.power_range_dBm    = exp.power_range;

    % ---- 初始化 ----
    rng(42);
    [SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(ov);
    [Bands, ~, Config] = init_m1_channel(Config, APs);
    [SpatialFP, ~] = init_m25_single_source_fp(APs, Bands, GridValid, Config, SourceTemplates);

    T = Config.m0.T_total;
    B = Config.m0.num_bands;
    M = APs.num;

    % ---- M0 循环 ----
    FrameStates = cell(1, T);
    for t = 1:T
        [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
            M0State, SourceTemplates, GridValid, Config, t, M0Logs);
        FrameStates{t} = FrameState_t;
    end

    % ---- M1 观测（覆写 n0） ----
    ChannelState = init_m1_channel_state(Config, APs, Bands);
    ChannelState.noise_n0_dBmHz = exp.n0_dBmHz * ones(1, B);

    Y_dBm_all = zeros(M, B, T);
    for t = 1:T
        [ChannelState, ObsFrame] = step_m1_generate_obs( ...
            FrameStates{t}, APs, Bands, ChannelState, Config, t);
        Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;
    end

    % ---- match_cfg（实验 C 启用约束 alpha） ----
    match_cfg = struct();
    match_cfg.mode = 'shifted_euclidean';
    if ei == 3
        % 库 ref_power = 0 dBm, target 范围 [-15, -5]
        % alpha = ref - tx, 所以 alpha in [5, 15]
        match_cfg.alpha_min_dB = 3;    % 留裕量
        match_cfg.alpha_max_dB = 20;
    else
        % 固定 0 dBm, alpha 应接近 0
        match_cfg.alpha_min_dB = -3;
        match_cfg.alpha_max_dB = 3;
    end

    % ---- 直接 WKNN 定位（仅 Band 1/2 的 ordinary_target） ----
    errors = [];
    bands  = [];
    for t = 1:T
        for bi = 1:numel(test_bands)
            b = test_bands(bi);
            apb = FrameStates{t}.active_per_band(b);
            if ~apb.has_source, continue; end
            if ~strcmp(apb.source_type, 'ordinary_target'), continue; end

            F_obs = Y_dBm_all(:, b, t);
            [dv, ~] = compute_wknn_distance_power_corrected(F_obs, SpatialFP.band(b), match_cfg);
            [ni, ~, nw] = select_knn_neighbors_m4(dv, K);
            ep = estimate_position_wknn_m4(SpatialFP.grid_xy, ni, nw);

            err = norm(ep - apb.true_pos_xy);
            errors(end+1) = err; %#ok<AGROW>
            bands(end+1)  = b;   %#ok<AGROW>
        end
    end

    % ---- 统计 ----
    fprintf('\n  总体 (n=%d): Mean=%.2f m, RMSE=%.2f m, Median=%.2f m\n', ...
        numel(errors), mean(errors), sqrt(mean(errors.^2)), median(errors));

    for bi = 1:numel(test_bands)
        b = test_bands(bi);
        mask = (bands == b);
        if any(mask)
            eb = errors(mask);
            fprintf('  Band %d (n=%d): Mean=%.2f m, RMSE=%.2f m, Median=%.2f m\n', ...
                b, sum(mask), mean(eb), sqrt(mean(eb.^2)), median(eb));
        end
    end

    all_results{ei} = struct('name', exp.name, 'errors', errors, 'bands', bands);
end

%% ========== 对比可视化 ==========

figure('Name', 'Sanity Check A/B/C', 'Position', [100 100 1200 500]);

% --- CDF 对比 ---
subplot(1, 2, 1);
colors = [0 0.5 1; 0 0.8 0.3; 0.9 0.2 0.2];
hold on;
leg = {};
for ei = 1:3
    e_sorted = sort(all_results{ei}.errors);
    cdf_y = (1:numel(e_sorted)) / numel(e_sorted);
    plot(e_sorted, cdf_y, '-', 'Color', colors(ei,:), 'LineWidth', 2);
    leg{ei} = all_results{ei}.name;
end
hold off;
xlabel('定位误差 (m)'); ylabel('CDF');
title('实验 A/B/C 误差 CDF 对比');
legend(leg, 'Location', 'southeast'); grid on;

% --- 箱线图对比 ---
subplot(1, 2, 2);
grp_data  = [];
grp_label = [];
for ei = 1:3
    e = all_results{ei}.errors;
    grp_data  = [grp_data, e];                              %#ok<AGROW>
    grp_label = [grp_label, ei * ones(1, numel(e))];        %#ok<AGROW>
end
boxplot(grp_data, grp_label, 'Labels', {'A', 'B', 'C'});
ylabel('定位误差 (m)');
title('实验 A/B/C 误差箱线图'); grid on;

sgtitle('Band 1/2 Sanity Check — 分阶段定位验证');

fprintf('\n====== Sanity Check 完成 ======\n');
