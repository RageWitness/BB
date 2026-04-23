%% CA_EVALUATE_CALIBRATION  标校前后指纹库定位效果对比
%
%  同 main_simulation 的场景/target 配置（相同 rng seed），
%  跑一次仿真 → M4(标校前) → CA 更新 → M4(标校后) → 对比分析。

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  CA 标校前后定位效果对比\n');
fprintf('============================================\n\n');

%% ========== 1. 初始化（与 main 相同） ==========

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
sim_override.fp.grid_step = 3;
sim_override.m4.fingerprint_type = 'centered_dBm';
sim_override.m4.fp_distance      = 'L2';
sim_override.debug.expose_true_source_state = true;

sim_override.m0.source.broadband_cal.schedule_mode = 'manual';
sim_override.m0.source.broadband_cal.manual_schedule = [
    struct('frame_range', [100, 150],  'pos_xy', [120, 150]);
    struct('frame_range', [300, 350],  'pos_xy', [75, 228]);
];

fp_cache_file = 'cache/SpatialFP_mwm_awgn.mat';

[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);
Config = CA_fill_defaults(Config);

T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;

%% ========== 2. 指纹库 ==========

fprintf('--- 构建/加载指纹库 ---\n');
if exist(fp_cache_file, 'file')
    S = load(fp_cache_file);
    SpatialFP = S.SpatialFP;
    fprintf('  从缓存加载: %s\n', fp_cache_file);
else
    [SpatialFP, ~] = init_m25_single_source_fp(APs, Bands, GridValid, Config, SourceTemplates);
end

%% ========== 3. 主循环 ==========

fprintf('--- 仿真主循环 (%d 帧) ---\n', T);
FrameStates = cell(1, T);
ObsFrames   = cell(1, T);
Y_dBm_all   = zeros(M, B, T);
Y_lin_all   = zeros(M, B, T);

tic;
for t = 1:T
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);
    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);
    FrameStates{t} = FrameState_t;
    ObsFrames{t}   = ObsFrame;
    Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;
    Y_lin_all(:, :, t) = ObsFrame.Y_lin;
end
fprintf('  主循环耗时: %.2f s\n', toc);

%% ========== 4. 构建 DSI → Labels → EventList ==========

fprintf('--- 构建 DSI + Labels + EventList ---\n');
DSI_all = cell(T, B);
for t = 1:T
    for b = 1:B
        DSI_all{t, b} = build_driven_source_input( ...
            FrameStates{t}.active_per_band(b), Config);
    end
end

ExternalLabelsRaw = build_labels_from_dsi(DSI_all, Config);
LabelTable = ingest_external_source_labels(ExternalLabelsRaw, Config);
SourceContext = bind_external_labels_to_events(LabelTable, FrameStates, Config);
EventList = local_build_event_list(SourceContext, Y_dBm_all, Y_lin_all, Config);
fprintf('  EventList: %d 个事件\n', numel(EventList));

%% ========== 5. M4 定位 — 标校前 ==========

fprintf('\n===== M4 定位 — 标校前 =====\n');
OnlineSpatialFP = CA_init_online_spatialfp(SpatialFP);
LocResults_before = run_m4_wknn_localization(EventList, OnlineSpatialFP, FrameStates, Config);

%% ========== 6. CA 更新 ==========

fprintf('\n===== CA 标校更新（按帧） =====\n');
[OnlineSpatialFP_new, CAResult] = CA_run_update_frames( ...
    OnlineSpatialFP, SourceContext, Y_lin_all, Y_dBm_all, Config);

fprintf('  CA status:          %s\n', CAResult.status);
fprintf('  n_samples_total:    %d\n', CAResult.n_samples_total);
fprintf('  n_samples_valid:    %d\n', CAResult.n_samples_valid);
fprintf('  n_eff_per_band:     %s\n', mat2str(CAResult.n_eff_per_band, 2));
if ~isempty(CAResult.update_log)
    fprintf('  n_updated_entries:  %d\n', CAResult.update_log.n_updated_entries);
    fprintf('  mean_abs_delta:     %.4f dB\n', CAResult.update_log.mean_abs_delta);
    fprintf('  max_abs_delta:      %.4f dB\n', CAResult.update_log.max_abs_delta);
end

%% ========== 7. M4 定位 — 标校后 ==========

fprintf('\n===== M4 定位 — 标校后 =====\n');
LocResults_after = run_m4_wknn_localization(EventList, OnlineSpatialFP_new, FrameStates, Config);

%% ========== 8. 误差提取 ==========

[err_before, bands_before, types_before] = extract_errors(LocResults_before);
[err_after,  bands_after,  types_after]  = extract_errors(LocResults_after);

fprintf('\n===== 总体统计 =====\n');
fprintf('  标校前:  n=%d  Mean=%.2fm  Median=%.2fm  RMSE=%.2fm\n', ...
    numel(err_before), mean(err_before), median(err_before), sqrt(mean(err_before.^2)));
fprintf('  标校后:  n=%d  Mean=%.2fm  Median=%.2fm  RMSE=%.2fm\n', ...
    numel(err_after), mean(err_after), median(err_after), sqrt(mean(err_after.^2)));
fprintf('  变化:    %.2fm (%.1f%%)\n', ...
    mean(err_after) - mean(err_before), ...
    (mean(err_after) - mean(err_before)) / mean(err_before) * 100);

%% ========== 9. 分频带统计表 ==========

unique_bands = unique(bands_before);
fprintf('\n--- 分频带定位误差对比 ---\n');
fprintf('  %-6s | %-6s | %-22s | %-22s\n', ...
    'Band', 'N', 'Before(Mean/Med/RMSE)', 'After(Mean/Med/RMSE)');
for bi = 1:numel(unique_bands)
    bb = unique_bands(bi);
    eb = err_before(bands_before == bb);
    ea = err_after(bands_after == bb);
    fprintf('  %-6d | %-6d | %5.1f / %5.1f / %5.1f | %5.1f / %5.1f / %5.1f\n', ...
        bb, numel(eb), ...
        mean(eb), median(eb), sqrt(mean(eb.^2)), ...
        mean(ea), median(ea), sqrt(mean(ea.^2)));
end

%% ========== 图1：场景布局 + 定位对比（前/后） ==========

figure('Name', 'CA 标校前后对比', 'Position', [50 50 1200 550]);

for pass = 1:2
    subplot(1, 2, pass);
    switch pass
        case 1, LR = LocResults_before; tag = '标校前';
        case 2, LR = LocResults_after;  tag = '标校后';
    end

    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    for k = 1:numel(LR)
        lr = LR(k);
        if isempty(lr.true_pos_xy) || all(lr.true_pos_xy == 0), continue; end
        switch lr.type_hat
            case 'target'
                clr = [0.9 0.2 0.2];
            case 'opportunistic'
                clr = [1.0 0.7 0.0];
            otherwise
                clr = [0.5 0.5 0.5];
        end
        plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', ...
            'Color', clr, 'MarkerSize', 8, 'MarkerFaceColor', clr);
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
            'Color', clr*0.6, 'MarkerSize', 10, 'LineWidth', 2);
        plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
             [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
             '-', 'Color', min(1, clr*0.65 + 0.25), 'LineWidth', 1);
    end

    xlabel('X (m)'); ylabel('Y (m)');
    title(tag);
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);

    leg_h = gobjects(5, 1);
    leg_h(1) = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    leg_h(2) = plot(NaN, NaN, 'o', 'Color', [0.9 0.2 0.2], 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.2 0.2]);
    leg_h(3) = plot(NaN, NaN, 'x', 'Color', [0.9 0.2 0.2]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    leg_h(4) = plot(NaN, NaN, 'o', 'Color', [1 0.7 0], 'MarkerSize', 8, 'MarkerFaceColor', [1 0.7 0]);
    leg_h(5) = plot(NaN, NaN, 'x', 'Color', [1 0.7 0]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    legend(leg_h, {'AP', 'target true', 'target est', 'opp true', 'opp est'}, ...
        'Location', 'best');
end
sgtitle('M4 WKNN 定位结果 — 标校前 / 标校后');

%% ========== 图2：分频带误差统计 (2x3 subplot) ==========

figure('Name', '分频带定位误差对比', 'Position', [50 50 1400 900]);
band_colors = [0.2 0.6 1; 0 0.8 0.4; 1 0.7 0; 0.9 0.2 0.2];
n_ub = numel(unique_bands);

% --- 2.1 箱线图 ---
subplot(2, 3, 1);
grp_data = []; grp_label = {};
for bi = 1:n_ub
    bb = unique_bands(bi);
    eb = err_before(bands_before == bb);
    ea = err_after(bands_after == bb);
    grp_data  = [grp_data, eb, ea]; %#ok<AGROW>
    grp_label = [grp_label, ...
        repmat({sprintf('B%d 前', bb)}, 1, numel(eb)), ...
        repmat({sprintf('B%d 后', bb)}, 1, numel(ea))]; %#ok<AGROW>
end
if numel(grp_data) > 1
    boxplot(grp_data, grp_label, 'GroupOrder', unique(grp_label, 'stable'));
end
ylabel('定位误差 (m)');
title('各频带误差箱线图');
grid on;

% --- 2.2 CDF 对比 ---
subplot(2, 3, 2);
hold on;
leg_items = {};
for bi = 1:n_ub
    bb = unique_bands(bi);
    ci = min(bb, size(band_colors, 1));

    eb = sort(err_before(bands_before == bb));
    ea = sort(err_after(bands_after == bb));

    if ~isempty(eb)
        plot(eb, (1:numel(eb))/numel(eb), '--', 'Color', band_colors(ci,:), 'LineWidth', 1.5);
        leg_items{end+1} = sprintf('B%d 前 (med=%.1f)', bb, median(eb)); %#ok<AGROW>
    end
    if ~isempty(ea)
        plot(ea, (1:numel(ea))/numel(ea), '-', 'Color', band_colors(ci,:), 'LineWidth', 2);
        leg_items{end+1} = sprintf('B%d 后 (med=%.1f)', bb, median(ea)); %#ok<AGROW>
    end
end
hold off;
xlabel('定位误差 (m)'); ylabel('CDF');
title('各频带误差 CDF (虚线=前, 实线=后)');
legend(leg_items, 'Location', 'southeast', 'FontSize', 6);
grid on;

% --- 2.3 Mean/RMSE 柱状图 ---
subplot(2, 3, 3);
bar_data = zeros(n_ub, 4);
for bi = 1:n_ub
    bb = unique_bands(bi);
    eb = err_before(bands_before == bb);
    ea = err_after(bands_after == bb);
    bar_data(bi, :) = [mean(eb), sqrt(mean(eb.^2)), ...
                        mean(ea), sqrt(mean(ea.^2))];
end
bh = bar(unique_bands, bar_data);
bh(1).FaceColor = [0.7 0.85 1];   bh(2).FaceColor = [1 0.7 0.7];
bh(3).FaceColor = [0.3 0.6 0.9];  bh(4).FaceColor = [0.9 0.3 0.3];
xlabel('频带'); ylabel('误差 (m)');
title('Mean / RMSE 对比');
legend('Mean(前)', 'RMSE(前)', 'Mean(后)', 'RMSE(后)', 'Location', 'best', 'FontSize', 6);
grid on;

% --- 2.4 按源类型分频带误差 ---
subplot(2, 3, 4);
type_list  = {'target', 'opportunistic'};
type_colors = [0.9 0.2 0.2; 1 0.7 0];
hold on;
leg_h2 = gobjects(0); leg_s2 = {};
for ti = 1:numel(type_list)
    tmask_b = strcmp(types_before, type_list{ti});
    tmask_a = strcmp(types_after,  type_list{ti});
    for bi = 1:n_ub
        bb = unique_bands(bi);
        bmask_b = bands_before == bb;
        bmask_a = bands_after == bb;

        eb = err_before(tmask_b & bmask_b);
        ea = err_after(tmask_a & bmask_a);

        x_b = bb + (ti - 1.5)*0.25 - 0.08;
        x_a = bb + (ti - 1.5)*0.25 + 0.08;

        if ~isempty(eb)
            plot(x_b * ones(size(eb)), eb, 'o', 'Color', type_colors(ti,:)*0.6, ...
                'MarkerSize', 4);
            plot(x_b + [-0.05 0.05], [mean(eb) mean(eb)], '-', ...
                'Color', type_colors(ti,:)*0.4, 'LineWidth', 2);
        end
        if ~isempty(ea)
            h = plot(x_a * ones(size(ea)), ea, 'o', 'Color', type_colors(ti,:), ...
                'MarkerSize', 5, 'MarkerFaceColor', type_colors(ti,:));
            plot(x_a + [-0.05 0.05], [mean(ea) mean(ea)], '-', ...
                'Color', type_colors(ti,:)*0.7, 'LineWidth', 2);
            if bi == 1
                leg_h2(end+1) = h; %#ok<AGROW>
                leg_s2{end+1} = sprintf('%s(后)', type_list{ti}); %#ok<AGROW>
            end
        end
    end
end
hold off;
xlabel('频带'); ylabel('定位误差 (m)');
title('按源类型分频带 (空心=前, 实心=后)');
set(gca, 'XTick', unique_bands);
if ~isempty(leg_h2)
    legend(leg_h2, leg_s2, 'Location', 'best');
end
grid on;

% --- 2.5 误差 vs 事件长度 ---
subplot(2, 3, 5);
loc_mask = arrayfun(@(ev) ...
    strcmp(ev.route_action,'localize_only') || ...
    strcmp(ev.route_action,'localize_then_calibrate'), EventList);
loc_ev = EventList(loc_mask);
valid_b_mask = ~arrayfun(@(r) isnan(r.loc_error), LocResults_before);
if numel(loc_ev) == sum(valid_b_mask)
    durations = [loc_ev.duration];
    scatter(durations, err_before, 25, bands_before, 'o');
    hold on;
    scatter(durations, err_after, 25, bands_after, 'filled');
    hold off;
    colorbar; clim([0.5 B+0.5]);
    xlabel('事件持续帧数'); ylabel('定位误差 (m)');
    title('误差 vs 事件长度 (空心=前, 实心=后)');
    grid on;
end

% --- 2.6 误差 vs 功率 ---
subplot(2, 3, 6);
if numel(loc_ev) == sum(valid_b_mask)
    powers = [loc_ev.power_level_est];
    scatter(powers, err_before, 25, bands_before, 'o');
    hold on;
    scatter(powers, err_after, 25, bands_after, 'filled');
    hold off;
    colorbar; clim([0.5 B+0.5]);
    xlabel('事件功率等级 (dBm)'); ylabel('定位误差 (m)');
    title('误差 vs 功率 (空心=前, 实心=后)');
    grid on;
end

sgtitle('M4 分频带定位误差 — 标校前/后 对比');

%% ========== 图3：指纹库变化热力图 ==========

figure('Name', '指纹库变化热力图', 'Position', [100 100 1200 400]);
for b = 1:min(B, 4)
    subplot(1, min(B,4), b);
    delta_map = OnlineSpatialFP_new.band(b).centered_dBm - SpatialFP.band(b).centered_dBm;
    mean_delta = mean(delta_map, 1);
    scatter(GridValid.xy(:,1), GridValid.xy(:,2), 12, mean_delta, 'filled');
    colorbar;
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('Band %d CA更新量 (dB)', b));
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);
end
sgtitle('CA 标校后指纹库变化 (越接近0表示未修改)');

%% ========== 图4：逐事件误差变化 ==========

figure('Name', '逐事件误差变化', 'Position', [100 550 900 400]);
n_ev = min(numel(err_before), numel(err_after));
delta_err = err_after(1:n_ev) - err_before(1:n_ev);
bar_colors = zeros(n_ev, 3);
bar_colors(delta_err < 0, :) = repmat([0.2 0.7 0.3], sum(delta_err < 0), 1);
bar_colors(delta_err >= 0, :) = repmat([0.9 0.3 0.3], sum(delta_err >= 0), 1);
for i = 1:n_ev
    bar(i, delta_err(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
    hold on;
end
hold off;
xlabel('事件序号'); ylabel('Δ误差 (m), 负=改善');
title(sprintf('逐事件误差变化 (改善: %d/%d, %.0f%%)', ...
    sum(delta_err < 0), n_ev, sum(delta_err < 0)/n_ev*100));
grid on;

%% ========== CA 诊断输出 ==========

if ~isempty(CAResult.delta)
    fprintf('\n--- Lambda 分布 ---\n');
    for b = 1:B
        lam = CAResult.delta.band(b).lambda;
        valid = CAResult.delta.band(b).valid;
        lam_v = lam(valid);
        if ~isempty(lam_v)
            fprintf('  band %d: mean=%.4f  std=%.4f  min=%.4f  max=%.4f\n', ...
                b, mean(lam_v), std(lam_v), min(lam_v), max(lam_v));
        end
    end
end

frob_sq = 0;
for b = 1:B
    d = OnlineSpatialFP_new.band(b).centered_dBm - SpatialFP.band(b).centered_dBm;
    frob_sq = frob_sq + sum(d(:).^2);
end
fprintf('\n||centered_after - centered_before||_F = %.4f\n', sqrt(frob_sq));

fprintf('\n============================================\n');
fprintf('  CA 评估完成\n');
fprintf('============================================\n');


%% ==================== 辅助函数 ====================

function [errs, bands, types] = extract_errors(LR)
    valid_mask = ~arrayfun(@(r) isnan(r.loc_error), LR);
    vlr = LR(valid_mask);
    errs  = [vlr.loc_error];
    bands = [vlr.band_id];
    types = {vlr.type_hat};
end


function EventList = local_build_event_list(SourceContext, Y_dBm_all, Y_lin_all, Config)
    LC = source_label_constants();
    EventList = [];
    if isempty(SourceContext.label_table), return; end
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
        ts = lbl.start_frame; te = lbl.end_frame;
        ev.obs_segment_dBm = squeeze(Y_dBm_all(:, b, ts:te));
        ev.obs_segment_lin = squeeze(Y_lin_all(:, b, ts:te));
        if ev.duration == 1
            ev.obs_segment_dBm = ev.obs_segment_dBm(:);
            ev.obs_segment_lin = ev.obs_segment_lin(:);
        end

        ev.label = lbl.label;
        if LC.name_map.isKey(lbl.label), ev.type_hat = LC.name_map(lbl.label);
        else, ev.type_hat = 'unknown'; end
        if LC.route_map.isKey(lbl.label), ev.route_action = LC.route_map(lbl.label);
        else, ev.route_action = 'hold'; end

        ev.linked_template_key = lbl.template_key;
        ev.source_uid = lbl.source_uid;
        ev.hold_reason = '';
        ev.upgrade_hint = 'none';

        ev.location_prior = lbl.location_prior;
        ev.power_prior    = lbl.power_prior;
        ev.metadata       = lbl.metadata;

        if isfield(ev.metadata, 'source_type_name')
            ev.source_type_name = ev.metadata.source_type_name;
        else
            ev.source_type_name = ev.type_hat;
        end
        ev.location_prior_type = ev.location_prior.type;
        ev.power_prior_type    = ev.power_prior.type;
        if isfield(ev.metadata, 'instance_id'), ev.instance_id = ev.metadata.instance_id;
        else, ev.instance_id = 0; end
        ev.is_calibration_source  = (lbl.label == LC.PERSISTENT_CAL || lbl.label == LC.BROADBAND_CAL);
        ev.is_target_source       = (lbl.label == LC.TARGET);
        ev.is_opportunistic_source = (lbl.label == LC.OPPORTUNISTIC);

        P_bar = mean(ev.obs_segment_dBm, 2);
        ev.power_level_est = mean(P_bar);
        if ev.duration > 1
            ev.power_stability_est = var(mean(ev.obs_segment_dBm, 1));
        else
            ev.power_stability_est = 0;
        end

        ev.score_trusted = 0; ev.score_prior_pos = 0;
        ev.score_prior_time = 0; ev.score_target = 0;
        ev.band_coverage_vec = zeros(1, Config.m0.num_bands);
        ev.band_coverage_vec(b) = 1;
        ev.n_valid_ap = M_ap;

        if isempty(EventList), EventList = ev;
        else, EventList(end+1) = ev; end %#ok<AGROW>
    end
end
