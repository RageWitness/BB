%% CA_ALGORITHMA_MANUAL_BROADBAND_DEMO
% Trusted-source RF_raw calibration demo with manual broadband sources.
%
% This script does not modify main_simulation.m. It builds a manual
% broadband-calibration scenario, runs localization with the original
% SpatialFP, builds a pending calibrated RF_raw map, and runs localization
% again with that pending map.

clear; clc; close all;
rng(42);

fprintf('============================================\n');
fprintf('  CA Algorithm A: manual broadband calibration demo\n');
fprintf('============================================\n\n');

%% ========== User entry ==========

% Manual broadband calibration source setting.
% Each row is one broadband calibration source segment.
manual_bc_pos_xy = [
    120, 150;
     75, 228
];
manual_bc_start_frame = [110, 300];
manual_bc_duration_frames = [40, 40];
manual_bc_tx_power_dBm = 20;

% Simulation and fingerprint cache.
fp_cache_enable = true;
fp_cache_file = fullfile('cache', 'SpatialFP_lognormal.mat');

% M4 setting. The calibration target is RF_raw.
m4_fingerprint_type = 'rf_raw';
m4_fp_distance = 'L2';

% CA Algorithm A parameter override.
ca_cfg = struct();
ca_cfg.N_min = 1;
ca_cfg.N_max = 5;
ca_cfg.K_coarse = 100;
ca_cfg.TolX = 1e-6;
ca_cfg.fusion_mode = 'mean';   % mean | median
ca_cfg.verbose = true;

%% ========== Build override ==========

sim_override = struct();

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
sim_override.m4.fingerprint_type = m4_fingerprint_type;
sim_override.m4.fp_distance = m4_fp_distance;
sim_override.debug.expose_true_source_state = true;

sim_override.m0.source.broadband_cal.schedule_mode = 'manual';
sim_override.m0.source.broadband_cal.manual_schedule = ...
    build_manual_schedule(manual_bc_pos_xy, manual_bc_start_frame, manual_bc_duration_frames);
sim_override.m0.source.broadband_cal.count = size(manual_bc_pos_xy, 1);
sim_override.m0.source.broadband_cal.tx_power_dBm = manual_bc_tx_power_dBm;

sim_override.ca.trusted_rss = ca_cfg;

%% ========== Init ==========

[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

if fp_cache_enable && exist(fp_cache_file, 'file')
    fprintf('[Cache] Load SpatialFP: %s\n', fp_cache_file);
    S = load(fp_cache_file);
    SpatialFP = S.SpatialFP;
    SignatureLib = S.SignatureLib; %#ok<NASGU>
else
    [SpatialFP, SignatureLib] = init_m25_single_source_fp( ... %#ok<ASGLU>
        APs, Bands, GridValid, Config, SourceTemplates);
    if fp_cache_enable
        cache_dir = fileparts(fp_cache_file);
        if ~isempty(cache_dir) && ~exist(cache_dir, 'dir')
            mkdir(cache_dir);
        end
        save(fp_cache_file, 'SpatialFP', 'SignatureLib', '-v7.3');
        fprintf('[Cache] Save SpatialFP: %s\n', fp_cache_file);
    end
end

T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;

fprintf('\n--- Simulation setting ---\n');
fprintf('  T=%d, B=%d, M=%d, Grid=%d\n', T, B, M, GridValid.Nvalid);
for k = 1:numel(sim_override.m0.source.broadband_cal.manual_schedule)
    ms = sim_override.m0.source.broadband_cal.manual_schedule(k);
    fprintf('  BC%d: frame=[%d,%d], pos=(%.1f, %.1f), power=%.1f dBm\n', ...
        k, ms.frame_range(1), ms.frame_range(2), ms.pos_xy(1), ms.pos_xy(2), ...
        manual_bc_tx_power_dBm);
end

%% ========== Run simulation ==========

FrameStates = cell(1, T);
ObsFrames = cell(1, T);
Y_dBm_all = zeros(M, B, T);
Y_lin_all = zeros(M, B, T);
occupancy_matrix = zeros(T, B);

fprintf('\n--- Main loop (%d frames) ---\n', T);
tic;
for t = 1:T
    [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
        M0State, SourceTemplates, GridValid, Config, t, M0Logs);

    [ChannelState, ObsFrame] = step_m1_generate_obs( ...
        FrameState_t, APs, Bands, ChannelState, Config, t);

    FrameStates{t} = FrameState_t;
    ObsFrames{t} = ObsFrame;
    Y_dBm_all(:, :, t) = ObsFrame.Y_dBm;
    Y_lin_all(:, :, t) = ObsFrame.Y_lin;

    for b = 1:B
        apb = FrameState_t.active_per_band(b);
        if apb.has_source
            switch apb.source_type
                case 'broadband_cal',  occupancy_matrix(t, b) = 1;
                case 'opportunistic',  occupancy_matrix(t, b) = 2;
                case 'persistent_cal', occupancy_matrix(t, b) = 3;
                case 'target',         occupancy_matrix(t, b) = 4;
            end
        end
    end
end
fprintf('Simulation done: %.2f s\n', toc);

%% ========== Labels and EventList ==========

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
EventList = local_build_event_list_from_context(SourceContext, Y_dBm_all, Y_lin_all, Config);

%% ========== Before localization ==========

fprintf('\n===== M4 with original SpatialFP =====\n');
LocResults_before = run_m4_wknn_localization(EventList, SpatialFP, FrameStates, Config);

%% ========== CA Algorithm A ==========

fprintf('\n===== CA Algorithm A trusted RSS calibration =====\n');
[PendingSpatialFP, CAResult] = CA_calibrate_spatialfp_trusted_sources( ...
    SpatialFP, SourceContext, Y_dBm_all, APs, Config);

fprintf('CA status: %s\n', CAResult.status);
fprintf('CA valid candidates: %d / %d\n', CAResult.n_candidates_valid, CAResult.n_candidates_total);
for b = 1:B
    fprintf('  Band %d: %s, candidates=%d, mean_delta=%.3f dB, n_hat=%s\n', ...
        b, CAResult.band(b).status, CAResult.band(b).n_candidates, ...
        CAResult.band(b).mean_delta_dB, mat2str(CAResult.band(b).n_hat_values, 3));
end

%% ========== After localization ==========

fprintf('\n===== M4 with pending calibrated SpatialFP =====\n');
LocResults_after = run_m4_wknn_localization(EventList, PendingSpatialFP, FrameStates, Config);

%% ========== Report and plots ==========

print_error_compare(LocResults_before, LocResults_after);
plot_occupancy_timeline(occupancy_matrix, B, T);
plot_scene_layout(GridValid, APs, M0Logs, Config);
plot_before_after_localization( ...
    LocResults_before, LocResults_after, GridValid, APs, Config);
plot_ca_delta_maps(SpatialFP, PendingSpatialFP, GridValid, Config);

fprintf('\n============================================\n');
fprintf('  CA Algorithm A demo done\n');
fprintf('============================================\n');


%% ==================== helpers ====================

function schedule = build_manual_schedule(pos_xy, start_frames, duration_frames)
    n = size(pos_xy, 1);
    if numel(start_frames) ~= n || numel(duration_frames) ~= n
        error('manual_bc_start_frame and manual_bc_duration_frames must match pos rows');
    end
    schedule = repmat(struct('frame_range', [1, 1], 'pos_xy', [0, 0]), n, 1);
    for i = 1:n
        sf = round(start_frames(i));
        ef = sf + round(duration_frames(i)) - 1;
        schedule(i).frame_range = [sf, ef];
        schedule(i).pos_xy = pos_xy(i, :);
    end
end


function EventList = local_build_event_list_from_context(SourceContext, Y_dBm_all, Y_lin_all, Config)
    LC = source_label_constants();
    EventList = [];
    if isempty(SourceContext.label_table)
        fprintf('[BuildEvents] empty label table\n');
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

        ev.location_prior = get_struct_field(lbl, 'location_prior', struct('type', 'none', 'value', []));
        ev.power_prior = get_struct_field(lbl, 'power_prior', struct('type', 'none', 'value', []));
        ev.metadata = get_struct_field(lbl, 'metadata', struct());

        if isfield(ev.metadata, 'source_type_name') && ~isempty(ev.metadata.source_type_name)
            ev.source_type_name = ev.metadata.source_type_name;
        else
            ev.source_type_name = ev.type_hat;
        end
        ev.location_prior_type = ev.location_prior.type;
        ev.power_prior_type = ev.power_prior.type;
        if isfield(ev.metadata, 'instance_id')
            ev.instance_id = ev.metadata.instance_id;
        else
            ev.instance_id = 0;
        end
        ev.is_calibration_source = (lbl.label == LC.PERSISTENT_CAL || lbl.label == LC.BROADBAND_CAL);
        ev.is_target_source = (lbl.label == LC.TARGET);
        ev.is_opportunistic_source = (lbl.label == LC.OPPORTUNISTIC);

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
    fprintf('[BuildEvents] %d labels -> %d events\n', ...
        numel(SourceContext.label_table), numel(EventList));
end


function v = get_struct_field(s, fname, default_val)
    if isfield(s, fname)
        v = s.(fname);
    else
        v = default_val;
    end
end


function print_error_compare(before, after)
    eb = collect_errors(before);
    ea = collect_errors(after);
    fprintf('\n--- Localization error compare ---\n');
    print_one_error_line('Before', eb);
    print_one_error_line('After ', ea);
    if ~isempty(eb) && ~isempty(ea)
        n = min(numel(eb), numel(ea));
        fprintf('  Delta mean: %.2f m\n', mean(ea(1:n)) - mean(eb(1:n)));
    end
end


function errs = collect_errors(LR)
    if isempty(LR)
        errs = [];
        return;
    end
    errs = [LR.loc_error];
    errs = errs(isfinite(errs));
end


function print_one_error_line(tag, errs)
    if isempty(errs)
        fprintf('  %s: no valid errors\n', tag);
        return;
    end
    fprintf('  %s: n=%d mean=%.2f med=%.2f rmse=%.2f p90=%.2f\n', ...
        tag, numel(errs), mean(errs), median(errs), sqrt(mean(errs.^2)), prctile(errs, 90));
end


function plot_occupancy_timeline(occupancy_matrix, B, T)
    type_names = {'空闲', 'broadband_cal', 'opportunistic', 'persistent_cal', 'target'};
    cmap = [1 1 1; 0.2 0.6 1; 0 0.8 0.4; 1 0.7 0; 0.9 0.2 0.2];
    figure('Name', 'CA-A 频带占用时间线', 'Position', [50 50 1200 350]);
    imagesc(1:T, 1:B, occupancy_matrix');
    colormap(cmap);
    cbar = colorbar;
    cbar.Ticks = [0, 1, 2, 3, 4];
    cbar.TickLabels = type_names;
    xlabel('帧号'); ylabel('频带');
    title('M0 频带占用时间线');
    set(gca, 'YTick', 1:B);
    grid on;
end


function plot_scene_layout(GridValid, APs, M0Logs, Config)
    figure('Name', 'CA-A 场景布局', 'Position', [700 50 700 600]);
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.85 0.85 0.85], 'MarkerSize', 3);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    h_bc = []; h_pc = []; h_opp = []; h_tg = [];
    if ~isempty(M0Logs.TruthLogAll)
        all_types = {M0Logs.TruthLogAll.source_type};
        all_xy = reshape([M0Logs.TruthLogAll.true_pos_xy], 2, [])';
        all_ids = [M0Logs.TruthLogAll.instance_id];
        h_bc = plot_unique_sources(all_xy, all_ids, strcmp(all_types, 'broadband_cal'), 'p', [0.2 0.5 1], 16);
        h_pc = plot_unique_sources(all_xy, all_ids, strcmp(all_types, 'persistent_cal'), 's', [0 0.8 0.4], 10);
        h_opp = plot_unique_sources(all_xy, all_ids, strcmp(all_types, 'opportunistic'), 'd', [1 0.7 0], 9);
        h_tg = plot_unique_sources(all_xy, all_ids, strcmp(all_types, 'target'), 'o', [0.1 0.2 0.2], 7);
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
    if ~isempty(h_bc), leg_h{end+1} = h_bc; leg_s{end+1} = 'broadband\_cal'; end
    if ~isempty(h_pc), leg_h{end+1} = h_pc; leg_s{end+1} = 'persistent\_cal'; end
    if ~isempty(h_opp), leg_h{end+1} = h_opp; leg_s{end+1} = 'opportunistic'; end
    if ~isempty(h_tg), leg_h{end+1} = h_tg; leg_s{end+1} = 'target'; end
    legend([leg_h{:}], leg_s, 'Location', 'best');
    hold off;
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


function plot_before_after_localization(before, after, GridValid, APs, Config)
    figure('Name', 'CA-A M4 定位结果对比', 'Position', [50 100 1450 650]);
    subplot(1,2,1);
    plot_loc_panel(before, GridValid, APs, Config, '标校前');
    subplot(1,2,2);
    plot_loc_panel(after, GridValid, APs, Config, '标校后');
    sgtitle('M4 WKNN 定位结果 — 原库 / CA-A 标校后库');
end


function plot_loc_panel(LR, GridValid, APs, Config, ttl)
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    for k = 1:numel(LR)
        lr = LR(k);
        if isempty(lr.true_pos_xy) || all(lr.true_pos_xy == 0)
            continue;
        end
        [clr, short_name] = color_for_type(lr.type_hat);
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
    title(ttl);
    axis equal; grid on;
    xlim(Config.area.x_range); ylim(Config.area.y_range);

    leg_h = gobjects(5, 1);
    leg_h(1) = plot(NaN, NaN, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    leg_h(2) = plot(NaN, NaN, 'o', 'Color', [0.9 0.2 0.2], 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.2 0.2]);
    leg_h(3) = plot(NaN, NaN, 'x', 'Color', [0.9 0.2 0.2]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    leg_h(4) = plot(NaN, NaN, 'o', 'Color', [1 0.7 0], 'MarkerSize', 8, 'MarkerFaceColor', [1 0.7 0]);
    leg_h(5) = plot(NaN, NaN, 'x', 'Color', [1 0.7 0]*0.6, 'MarkerSize', 10, 'LineWidth', 2);
    legend(leg_h, {'AP', 'target true', 'target est', 'opportunistic true', 'opportunistic est'}, ...
        'Location', 'best');
    hold off;
end


function [clr, short_name] = color_for_type(type_hat)
    switch type_hat
        case 'target'
            clr = [0.9 0.2 0.2];
            short_name = 'target';
        case 'opportunistic'
            clr = [1.0 0.7 0.0];
            short_name = 'opp';
        otherwise
            clr = [0.5 0.5 0.5];
            short_name = type_hat;
    end
end


function plot_ca_delta_maps(SpatialFP, PendingSpatialFP, GridValid, Config)
    B = SpatialFP.B;
    n_plot = min(B, 4);
    all_vals = [];
    mean_delta = cell(1, n_plot);
    for b = 1:n_plot
        D = get_rf_raw(PendingSpatialFP.band(b)) - get_rf_raw(SpatialFP.band(b));
        mean_delta{b} = mean(D, 1);
        all_vals = [all_vals; mean_delta{b}(:)]; %#ok<AGROW>
    end
    all_vals = all_vals(isfinite(all_vals));
    if isempty(all_vals)
        clim = [-1, 1];
    else
        a = max(abs(all_vals));
        if a <= 0, a = 1; end
        clim = [-a, a];
    end

    figure('Name', 'CA-A RF_raw 指纹库变化', 'Position', [80 150 1350 360]);
    for b = 1:n_plot
        subplot(1, n_plot, b);
        scatter(GridValid.xy(:,1), GridValid.xy(:,2), 12, mean_delta{b}, 'filled');
        caxis(clim); colorbar;
        xlabel('X (m)'); ylabel('Y (m)');
        title(sprintf('Band %d mean AP delta (dB)', b));
        axis equal; grid on;
        xlim(Config.area.x_range); ylim(Config.area.y_range);
    end
    sgtitle('CA-A 标校后 RF\_raw 指纹库变化：PendingSpatialFP - SpatialFP');
end


function F = get_rf_raw(band_fp)
    if isfield(band_fp, 'RF_raw')
        F = band_fp.RF_raw;
    else
        F = band_fp.F_dBm;
    end
end
