%% DEBUG_LOCALIZATION  诊断定位误差大的根因
%  在 main_simulation.m 运行完毕后执行（需要工作区中的变量）

fprintf('\n===== 定位误差诊断 =====\n\n');

%% 1. 检查 M3 事件检测：是否把多个源合并成了一个事件
fprintf('--- 1. 事件检测 vs 真值对比 ---\n');
for e = 1:min(numel(EventList), 10)
    ev = EventList(e);
    if ~strcmp(ev.route_action, 'localize_only') && ...
       ~strcmp(ev.route_action, 'localize_then_calibrate')
        continue;
    end

    b = ev.band_id;
    fprintf('\n  Event %d: band=%d, t=[%d,%d], type_hat=%s\n', ...
        ev.event_id, b, ev.t_start, ev.t_end, ev.type_hat);

    % 检查这个事件窗口内真值源的情况
    unique_positions = [];
    unique_types = {};
    for t = ev.t_start : ev.t_end
        fs = FrameStates{t};
        apb = fs.active_per_band(b);
        if apb.has_source
            pos = apb.true_pos_xy;
            stype = apb.source_type;
            fprintf('    t=%d: src=%s, pos=(%.1f,%.1f), power=%.1f dBm\n', ...
                t, stype, pos(1), pos(2), apb.tx_power_dBm);
            unique_positions(end+1,:) = pos; 
            unique_types{end+1} = stype; 
        else
            fprintf('    t=%d: 无源！\n', t);
        end
    end

    % 检查是否混入了多个不同位置的源
    if size(unique_positions, 1) > 1
        pos_spread = max(range(unique_positions, 1));
        n_unique_types = numel(unique(unique_types));
        fprintf('  >>> 位置离散度: %.1f m, 不同源类型数: %d\n', pos_spread, n_unique_types);
    end
end

%% 2. 逐事件验证 WKNN：拿真值位置直接查指纹库距离
fprintf('\n\n--- 2. WKNN 匹配诊断 ---\n');
for k = 1:min(numel(LocResults), 10)
    lr = LocResults(k);
    if isnan(lr.loc_error), continue; end

    b = lr.band_id;
    true_pos = lr.true_pos_xy;
    est_pos  = lr.est_pos_xy;

    % 真值位置最近的网格点
    d_to_grid = sqrt(sum((SpatialFP.grid_xy - true_pos).^2, 2));
    [d_nearest_grid, g_true] = min(d_to_grid);

    % 该网格点在指纹库中的中心化指纹
    F_db_centered_true = SpatialFP.band(b).centered_dBm(:, g_true);

    % 事件聚合指纹
    F_obs = lr.obs_fp_dBm;
    mu_obs = mean(F_obs);
    F_obs_centered = F_obs - mu_obs;

    % 真值网格点的功率修正距离
    d_true_grid = norm(F_obs_centered - F_db_centered_true);

    % 最佳匹配网格点
    g_best = lr.neighbor_idx(1);
    d_best = lr.neighbor_dist(1);
    best_pos = SpatialFP.grid_xy(g_best, :);

    fprintf('\n  Event %d (band %d):\n', lr.event_id, b);
    fprintf('    真值位置:     (%.1f, %.1f)\n', true_pos(1), true_pos(2));
    fprintf('    最近网格点:   g=%d, 距真值 %.2f m\n', g_true, d_nearest_grid);
    fprintf('    真值网格距离: d=%.4f (指纹空间)\n', d_true_grid);
    fprintf('    最佳匹配:     g=%d (%.1f, %.1f), d=%.4f\n', ...
        g_best, best_pos(1), best_pos(2), d_best);
    fprintf('    估计位置:     (%.1f, %.1f)\n', est_pos(1), est_pos(2));
    fprintf('    定位误差:     %.2f m\n', lr.loc_error);

    % 检查观测指纹与指纹库的数值范围
    F_db_raw_true = SpatialFP.band(b).F_dBm(:, g_true);
    fprintf('    观测指纹均值: %.1f dBm\n', mu_obs);
    fprintf('    指纹库均值:   %.1f dBm (g_true)\n', SpatialFP.band(b).mean_dBm(g_true));
    fprintf('    观测指纹范围: [%.1f, %.1f] dBm\n', min(F_obs), max(F_obs));
    fprintf('    指纹库范围:   [%.1f, %.1f] dBm (g_true)\n', min(F_db_raw_true), max(F_db_raw_true));

    % 中心化后对比
    fprintf('    |F_obs_c|:    %.4f\n', norm(F_obs_centered));
    fprintf('    |F_db_c|:     %.4f\n', norm(F_db_centered_true));
    fprintf('    cos相似度:    %.6f\n', ...
        dot(F_obs_centered, F_db_centered_true) / ...
        (norm(F_obs_centered) * norm(F_db_centered_true) + 1e-20));

    % 排名：真值网格点在距离排序中排第几
    ev_idx = find([EventList.event_id] == lr.event_id);
    ev = EventList(ev_idx);
    F_obs_recomp = aggregate_event_fingerprint_m4(ev);
    [dist_all, ~] = compute_wknn_distance_power_corrected(F_obs_recomp, SpatialFP.band(b));
    rank_true = sum(dist_all <= dist_all(g_true));
    fprintf('    真值网格排名: %d / %d\n', rank_true, SpatialFP.G);
end

%% 3. 单帧直接定位测试（跳过 M3，直接用真值帧做 WKNN）
fprintf('\n\n--- 3. 跳过M3，直接用真值帧定位 ---\n');
K = 5;
direct_errors = [];
direct_bands  = [];
for t = 1:10:T  % 每 10 帧取一个
    for b = 1:B
        apb = FrameStates{t}.active_per_band(b);
        if ~apb.has_source, continue; end
        if ~strcmp(apb.source_type, 'ordinary_target'), continue; end

        % 直接取这一帧的观测
        F_obs_direct = Y_dBm_all(:, b, t);

        % WKNN 定位
        [dv, ~] = compute_wknn_distance_power_corrected(F_obs_direct, SpatialFP.band(b));
        [ni, nd, nw] = select_knn_neighbors_m4(dv, K);
        ep = estimate_position_wknn_m4(SpatialFP.grid_xy, ni, nw);

        err = norm(ep - apb.true_pos_xy);
        direct_errors(end+1) = err; %#ok<AGROW>
        direct_bands(end+1)  = b;   %#ok<AGROW>

        if numel(direct_errors) <= 5
            fprintf('  t=%d, b=%d: true=(%.1f,%.1f), est=(%.1f,%.1f), err=%.2f m\n', ...
                t, b, apb.true_pos_xy(1), apb.true_pos_xy(2), ep(1), ep(2), err);
        end
    end
end

if ~isempty(direct_errors)
    fprintf('\n  直接定位统计 (n=%d):\n', numel(direct_errors));
    fprintf('    Mean: %.2f m\n', mean(direct_errors));
    fprintf('    RMSE: %.2f m\n', sqrt(mean(direct_errors.^2)));
    fprintf('    Median: %.2f m\n', median(direct_errors));

    % --- 分频带误差统计表 ---
    unique_bands = unique(direct_bands);
    fprintf('\n  分频带误差统计:\n');
    fprintf('  %-6s %-8s %-10s %-10s %-10s %-10s %-10s\n', ...
        'Band', 'N', 'Mean(m)', 'Median(m)', 'RMSE(m)', 'Max(m)', 'Min(m)');
    band_rmse = zeros(1, max(unique_bands));
    for ib = 1:numel(unique_bands)
        bid = unique_bands(ib);
        e_b = direct_errors(direct_bands == bid);
        band_rmse(bid) = sqrt(mean(e_b.^2));
        fprintf('  %-6d %-8d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n', ...
            bid, numel(e_b), mean(e_b), median(e_b), band_rmse(bid), max(e_b), min(e_b));
    end

    % --- 分频带误差图 (2x2) ---
    band_colors = lines(numel(unique_bands));
    figure('Name', '分频带定位误差', 'Position', [150 150 1100 800]);

    % 子图1：箱线图
    subplot(2,2,1);
    grp_data = [];
    grp_label = [];
    for ib = 1:numel(unique_bands)
        bid = unique_bands(ib);
        e_b = direct_errors(direct_bands == bid);
        grp_data  = [grp_data, e_b];          %#ok<AGROW>
        grp_label = [grp_label, bid * ones(1, numel(e_b))]; %#ok<AGROW>
    end
    boxplot(grp_data, grp_label);
    xlabel('Band'); ylabel('误差 (m)');
    title('各频带误差箱线图'); grid on;

    % 子图2：RMSE 柱状图
    subplot(2,2,2);
    bar_vals = band_rmse(unique_bands);
    bh = bar(1:numel(unique_bands), bar_vals, 0.5);
    bh.FaceColor = 'flat';
    for ib = 1:numel(unique_bands)
        bh.CData(ib,:) = band_colors(ib,:);
    end
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('Band %d', x), unique_bands, 'Uni', false));
    ylabel('RMSE (m)');
    title('各频带 RMSE'); grid on;

    % 子图3：CDF
    subplot(2,2,3);
    hold on;
    leg_str = cell(1, numel(unique_bands));
    for ib = 1:numel(unique_bands)
        bid = unique_bands(ib);
        e_b = sort(direct_errors(direct_bands == bid));
        cdf_y = (1:numel(e_b)) / numel(e_b);
        plot(e_b, cdf_y, '-', 'Color', band_colors(ib,:), 'LineWidth', 1.5);
        leg_str{ib} = sprintf('Band %d', bid);
    end
    hold off;
    xlabel('误差 (m)'); ylabel('CDF');
    title('各频带误差 CDF'); legend(leg_str, 'Location', 'southeast');
    grid on;

    % 子图4：散点图（按 Band 着色）
    subplot(2,2,4);
    hold on;
    for ib = 1:numel(unique_bands)
        bid = unique_bands(ib);
        idx_b = find(direct_bands == bid);
        scatter(idx_b, direct_errors(idx_b), 20, band_colors(ib,:), 'filled', 'MarkerFaceAlpha', 0.6);
    end
    hold off;
    xlabel('样本序号'); ylabel('误差 (m)');
    title('各频带误差散点'); legend(leg_str, 'Location', 'best');
    grid on;

    sgtitle('跳过M3 — 分频带直接定位误差分析');
end

%% 4. SNR 诊断
fprintf('\n\n--- 4. SNR 诊断 ---\n');

% 4a. 各频带噪声功率
fprintf('\n  [频带噪声功率]\n');
fprintf('  %-8s %-12s %-14s %-14s %-12s\n', ...
    'Band', 'Model', 'BW', 'n0(dBm/Hz)', 'N_floor(dBm)');
for b = 1:B
    bw = Bands.bw_Hz(b);
    if isfield(Config.m1, 'noise') && isfield(Config.m1.noise, 'n0_dBmHz')
        n0 = Config.m1.noise.n0_dBmHz(b);
    else
        n0 = -174;
    end
    N_dBm = n0 + 10 * log10(bw);
    fprintf('  %-8d %-12s %-14s %-14.1f %-12.1f\n', ...
        b, Bands.model{b}, format_bw(bw), n0, N_dBm);
end

% 4b. 逐帧逐频带采集 SNR 样本（用 ObsFrames.meta 中的真实接收信号功率）
fprintf('\n  [逐频带 SNR 统计 — 仅有源帧]\n');
fprintf('  SNR = rx_signal_dBm (接收信号功率,含衰落) - N_dBm (噪声总功率)\n\n');
snr_samples = cell(1, B);
rx_sig_samples = cell(1, B);
noise_samples = cell(1, B);
pl_samples = cell(1, B);

for t = 1:T
    for b = 1:B
        meta_b = ObsFrames{t}.meta.active_per_band(b);
        if ~meta_b.has_source, continue; end

        % 从 M1 记录的真实值取接收信号功率（含路径损耗+阴影，不含 AWGN）
        rx_sig = meta_b.rx_signal_dBm_per_ap;   % M x 1
        N_dBm  = meta_b.noise_power_dBm;        % 标量

        % 正确的 SNR：接收信号功率 / 噪声功率
        snr_per_ap = rx_sig(:) - N_dBm;          % M x 1, dB
        snr_samples{b}    = [snr_samples{b};    snr_per_ap];
        rx_sig_samples{b} = [rx_sig_samples{b}; rx_sig(:)];
        noise_samples{b}  = [noise_samples{b};  N_dBm];
        pl_samples{b}     = [pl_samples{b};     meta_b.pathloss_dB_per_ap(:)];
    end
end

fprintf('  %-6s %-10s %-12s %-10s %-10s %-10s %-14s %-14s %-14s\n', ...
    'Band', 'Model', 'N_noise(dBm)', 'N_samples', 'SNR_mean', 'SNR_med', ...
    'SNR<0(%)', 'RxSig_mean', 'PL_mean(dB)');
for b = 1:B
    if isempty(snr_samples{b}), continue; end
    snr = snr_samples{b};
    rx  = rx_sig_samples{b};
    pl  = pl_samples{b};
    % 噪声功率（取第一个值即可，同频带内固定）
    N_val = noise_samples{b}(1);
    fprintf('  %-6d %-10s %-12.1f %-10d %-10.1f %-10.1f %-14.1f %-14.1f %-14.1f\n', ...
        b, Bands.model{b}, N_val, numel(snr), mean(snr), median(snr), ...
        100*sum(snr < 0)/numel(snr), mean(rx), mean(pl));
end

% 4c. SNR 分布直方图
figure('Name', 'SNR 诊断', 'Position', [100 100 1200 500]);
for b = 1:B
    if isempty(snr_samples{b}), continue; end
    subplot(1, B, b);
    histogram(snr_samples{b}, 30, 'FaceColor', [0.3 0.6 0.9]);
    hold on;
    xline(0, 'r--', 'SNR=0', 'LineWidth', 2);
    hold off;
    xlabel('SNR (dB)'); ylabel('计数');
    title(sprintf('Band %d (%s)', b, Bands.model{b}));
    grid on;
end
sgtitle('各频带 SNR 分布 (有源帧, 全 AP)');

% 4d. 信号/噪声/观测 逐AP 示例（取第一个有源帧）
figure('Name', '信号-噪声-观测 分解', 'Position', [100 600 1400 400]);
for b = 1:B
    % 找第一个有源帧
    t_demo = 0;
    for t = 1:T
        if ObsFrames{t}.meta.active_per_band(b).has_source
            t_demo = t; break;
        end
    end
    if t_demo == 0, continue; end

    meta_b = ObsFrames{t_demo}.meta.active_per_band(b);
    rx_sig = meta_b.rx_signal_dBm_per_ap(:);  % 接收信号功率（含衰落，不含 AWGN）
    N_dBm  = meta_b.noise_power_dBm;
    Y_obs  = Y_dBm_all(:, b, t_demo);
    Y_obs  = Y_obs(:);

    % 逐 AP 的 SNR
    snr_ap = rx_sig - N_dBm;

    subplot(1, B, b);
    bar_data = [rx_sig, Y_obs];
    bh = bar(bar_data);
    bh(1).FaceColor = [0.3 0.7 0.3];
    bh(2).FaceColor = [0.9 0.3 0.3];
    hold on;
    yline(N_dBm, 'b--', sprintf('N=%.0fdBm', N_dBm), 'LineWidth', 1.5);
    hold off;
    xlabel('AP'); ylabel('dBm');
    title(sprintf('Band %d  P_{tx}=%.0fdBm  %s\nmean SNR=%.1f dB', ...
        b, meta_b.tx_power_dBm, Bands.model{b}, mean(snr_ap)));
    legend({'rx信号(含衰落)', '观测(信号+AWGN)', '噪声底'}, ...
        'Location', 'best', 'FontSize', 7);
    grid on;
    set(gca, 'XTick', 1:4:M);
end
sgtitle('各频带 接收信号 / AWGN观测 / 噪声底  逐AP分解');

fprintf('\n===== 诊断完成 =====\n');


%% ==================== 辅助函数 ====================

function s = format_bw(bw_hz)
    if bw_hz >= 1e6
        s = sprintf('%.0f MHz', bw_hz / 1e6);
    else
        s = sprintf('%.0f kHz', bw_hz / 1e3);
    end
end
