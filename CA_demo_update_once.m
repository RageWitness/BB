%% CA_DEMO_UPDATE_ONCE  独立 demo：跑一次完整 CA 更新
%
%  前置：运行过 main_simulation 生成 EventList 和 SpatialFP

clear; clc;
rng(42);
fprintf('========== CA_DEMO_UPDATE_ONCE ==========\n\n');

% --- 仿真配置 ---
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
sim_override.m0.T_total = 200;

% --- CA 参数覆盖（可选，演示用） ---
sim_override.ca.kmeans.K = 4;
sim_override.ca.gpr.learn_hyperparams = true;
sim_override.ca.kmeans.enable_blend = true;
sim_override.ca.update.enable_adaptive_lambda = true;

fprintf('[Demo] 初始化 M0 ...\n');
[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);
Config = CA_fill_defaults(Config);

T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;

% --- 跑仿真 ---
fprintf('[Demo] 跑 %d 帧仿真 ...\n', T);
FrameStates = cell(1, T);
ObsFrames   = cell(1, T);
Y_dBm_all   = zeros(M, B, T);
Y_lin_all   = zeros(M, B, T);

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

% --- DSI → Labels → EventList ---
fprintf('[Demo] 构建 DSI + Labels + EventList ...\n');
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

fprintf('[Demo] EventList: %d 个事件\n', numel(EventList));

% --- 构建指纹库 ---
fprintf('[Demo] 构建指纹库 ...\n');
if exist('cache/SpatialFP_mwm_Nl.mat', 'file')
    S = load('cache/SpatialFP_mwm_Nl.mat');
    SpatialFP = S.SpatialFP;
    fprintf('[Demo] 从缓存加载 SpatialFP\n');
else
    [SpatialFP, ~] = init_m25_single_source_fp(APs, Bands, GridValid, Config, SourceTemplates);
end

% --- 初始化在线地图 ---
OfflineSpatialFP = SpatialFP;
OnlineSpatialFP  = CA_init_online_spatialfp(OfflineSpatialFP);

% 保存 before 快照
before_snap = cell(1, B);
for b = 1:B
    before_snap{b} = OnlineSpatialFP.band(b).centered_dBm;
end

% --- 运行 CA 更新（按帧） ---
fprintf('\n========== 运行 CA 更新（按帧） ==========\n');
[OnlineSpatialFP_new, CAResult] = CA_run_update_frames( ...
    OnlineSpatialFP, SourceContext, Y_lin_all, Y_dBm_all, Config);

% --- 输出 ---
fprintf('\n========== CA 结果 ==========\n');
fprintf('  status:            %s\n', CAResult.status);
fprintf('  n_samples_total:   %d\n', CAResult.n_samples_total);
fprintf('  n_samples_valid:   %d\n', CAResult.n_samples_valid);

fprintf('\n--- per-band 帧样本统计 ---\n');
for b = 1:B
    fprintf('  band %d: n_eff=%.1f\n', b, CAResult.n_eff_per_band(b));
end

if ~isempty(CAResult.update_log)
    ul = CAResult.update_log;
    fprintf('\n  n_updated_entries: %d\n', ul.n_updated_entries);
    fprintf('  mean_abs_delta:    %.6f\n', ul.mean_abs_delta);
    fprintf('  max_abs_delta:     %.6f\n', ul.max_abs_delta);
end

% Frobenius norm
frob_sq = 0;
for b = 1:B
    diff_b = OnlineSpatialFP_new.band(b).centered_dBm - before_snap{b};
    frob_sq = frob_sq + sum(diff_b(:).^2);
end
frob_norm = sqrt(frob_sq);
fprintf('  ||centered_after - centered_before||_F: %.6f\n', frob_norm);

% OfflineSpatialFP 未变
offline_changed = false;
for b = 1:B
    if any(OfflineSpatialFP.band(b).centered_dBm(:) ~= SpatialFP.band(b).centered_dBm(:))
        offline_changed = true;
        break;
    end
end
fprintf('  OfflineSpatialFP 未改变: %s\n', mat2str(~offline_changed));
fprintf('  OnlineSpatialFP.update_round:          %d\n', OnlineSpatialFP_new.update_round);
fprintf('  numel(update_history):                 %d\n', numel(OnlineSpatialFP_new.update_history));

% --- lambda 统计 ---
if ~isempty(CAResult.delta)
    fprintf('\n--- Lambda 分布 ---\n');
    for b = 1:B
        lam = CAResult.delta.band(b).lambda;
        valid = CAResult.delta.band(b).valid;
        lam_v = lam(valid);
        if ~isempty(lam_v)
            fprintf('  band %d: mean=%.4f  std=%.4f  min=%.4f  max=%.4f\n', ...
                b, mean(lam_v), std(lam_v), min(lam_v), max(lam_v));
        else
            fprintf('  band %d: no valid entries\n', b);
        end
    end

    fprintf('\n--- per-band n_updated_entries ---\n');
    for b = 1:B
        n_upd_b = sum(CAResult.delta.band(b).valid(:));
        fprintf('  band %d: n_updated=%d\n', b, n_upd_b);
    end
end
if ~isempty(CAResult.models) && isfield(CAResult.models, 'band')
    fprintf('\n--- 模型诊断 ---\n');
    for b = 1:min(B, numel(CAResult.models.band))
        for m = 1:min(M, numel(CAResult.models.band(b).ap))
            for k = 1:Config.ca.kmeans.K
                mdl = CAResult.models.band(b).ap(m).cluster(k);
                fprintf('  (b=%d,m=%d,k=%d): n=%d  valid=%d  status=%s', ...
                    b, m, k, mdl.n_train, mdl.valid, mdl.status);
                if mdl.valid && ~isempty(fieldnames(mdl.learned_hyperparams))
                    hp = mdl.learned_hyperparams;
                    fprintf('  ell=[%.1f,%.1f] sf=%.2f sn=%.2f', ...
                        hp.ell_x, hp.ell_y, hp.sigma_f, hp.sigma_n);
                end
                fprintf('\n');
            end
        end
    end
end

fprintf('\n========== DEMO 完成 ==========\n');


%% ==================== 辅助函数 ====================

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
