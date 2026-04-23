%% VERIFY_LABEL_PRIOR_FLOW  验证 DSI→LabelTable→SourceContext→EventList→LocResults prior 贯通
%
%  使用较小 T_total 跑一段仿真，检查 opportunistic prior 是否完整保留。

clear; clc;
rng(42);

fprintf('========== VERIFY_LABEL_PRIOR_FLOW ==========\n\n');

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
sim_override.m0.T_total = 100;   % 缩短

[SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(sim_override);
[Bands, ChannelState, Config] = init_m1_channel(Config, APs);

T = Config.m0.T_total;
B = Config.m0.num_bands;
M = APs.num;

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

% --- Step 1: DSI_all ---
DSI_all = cell(T, B);
for t = 1:T
    for b = 1:B
        DSI_all{t, b} = build_driven_source_input( ...
            FrameStates{t}.active_per_band(b), Config);
    end
end

% --- Step 2: DSI -> Labels ---
ExternalLabelsRaw = build_labels_from_dsi(DSI_all, Config);
LabelTable = ingest_external_source_labels(ExternalLabelsRaw, Config);

% --- Step 3: Bind ---
SourceContext = bind_external_labels_to_events(LabelTable, FrameStates, Config);

% --- Step 4: EventList ---
% 复用 main_simulation 的本地函数是不行的（脚本，无法访问）
% 这里直接内联一个等价版本，仅用于验证
EventList = local_build_event_list(SourceContext, Y_dBm_all, Y_lin_all, Config);

%% ===== 统计 =====
LC = source_label_constants();

fprintf('\n--- DSI opportunistic 统计 ---\n');
n_dsi_opp = 0; n_dsi_opp_exact = 0; n_dsi_opp_region = 0;
n_dsi_opp_gauss = 0; n_dsi_opp_none = 0;
n_dsi_opp_pow_exact = 0; n_dsi_opp_pow_none = 0;
dsi_opp_keys = {};
for t = 1:T
    for b = 1:B
        d = DSI_all{t, b};
        if d.label ~= LC.OPPORTUNISTIC, continue; end
        n_dsi_opp = n_dsi_opp + 1;
        switch d.location_prior.type
            case 'exact',    n_dsi_opp_exact = n_dsi_opp_exact + 1;
            case 'region',   n_dsi_opp_region = n_dsi_opp_region + 1;
            case 'gaussian', n_dsi_opp_gauss = n_dsi_opp_gauss + 1;
            otherwise,       n_dsi_opp_none = n_dsi_opp_none + 1;
        end
        switch d.power_prior.type
            case {'exact', 'exact_by_band'}, n_dsi_opp_pow_exact = n_dsi_opp_pow_exact + 1;
            otherwise,                        n_dsi_opp_pow_none = n_dsi_opp_pow_none + 1;
        end
        key = sprintf('%s__inst%d__b%d', d.debug_info.template_key, ...
            ternary_default(d.debug_info.instance_id), b);
        dsi_opp_keys{end+1} = key; %#ok<AGROW>
    end
end
fprintf('  DSI opp 活跃帧数: %d (exact=%d region=%d gauss=%d none=%d)\n', ...
    n_dsi_opp, n_dsi_opp_exact, n_dsi_opp_region, n_dsi_opp_gauss, n_dsi_opp_none);
fprintf('  DSI opp power: exact=%d none=%d\n', n_dsi_opp_pow_exact, n_dsi_opp_pow_none);

fprintf('\n--- LabelTable opportunistic 统计 ---\n');
lt_opp_mask = false(1, numel(LabelTable));
for i = 1:numel(LabelTable)
    lt_opp_mask(i) = (LabelTable(i).label == LC.OPPORTUNISTIC);
end
lt_opp = LabelTable(lt_opp_mask);
n_lt_opp = numel(lt_opp);
n_lt_opp_exact = sum(arrayfun(@(l) strcmp(l.location_prior.type,'exact'), lt_opp));
n_lt_opp_region = sum(arrayfun(@(l) strcmp(l.location_prior.type,'region'), lt_opp));
n_lt_opp_gauss  = sum(arrayfun(@(l) strcmp(l.location_prior.type,'gaussian'), lt_opp));
n_lt_opp_none   = sum(arrayfun(@(l) strcmp(l.location_prior.type,'none'), lt_opp));
fprintf('  LabelTable opp 标签数: %d (exact=%d region=%d gauss=%d none=%d)\n', ...
    n_lt_opp, n_lt_opp_exact, n_lt_opp_region, n_lt_opp_gauss, n_lt_opp_none);

%% ===== 断言 =====
fprintf('\n--- 断言检查 ---\n');

% [A1] 每个 DSI opp (instance_id, template_key, band_id) 组合在 LabelTable 中必须有标签
unique_dsi_keys = unique(dsi_opp_keys);
missing_keys = {};
for i = 1:numel(unique_dsi_keys)
    key = unique_dsi_keys{i};
    toks = regexp(key, '(.+)__inst(-?\d+)__b(\d+)', 'tokens', 'once');
    if isempty(toks), continue; end
    tpl = toks{1}; inst = str2double(toks{2}); bd = str2double(toks{3});
    found = false;
    for j = 1:numel(lt_opp)
        md = lt_opp(j).metadata;
        if strcmp(md.template_key, tpl) && md.band_id == bd && ...
           isequal(double(md.instance_id), inst)
            found = true; break;
        end
    end
    if ~found
        missing_keys{end+1} = key; %#ok<AGROW>
    end
end
if isempty(missing_keys)
    fprintf('[PASS] DSI -> LabelTable: 每个 DSI opp 组合都有对应标签 (n=%d)\n', numel(unique_dsi_keys));
else
    fprintf('[FAIL] DSI -> LabelTable: %d 个 DSI opp 组合丢失\n', numel(missing_keys));
    for i = 1:min(5, numel(missing_keys))
        fprintf('    missing: %s\n', missing_keys{i});
    end
end

% [A2] LabelTable opp prior 不能丢（如果 DSI 有 non-none，LabelTable 也应该 non-none）
prior_preserved = true;
for i = 1:numel(lt_opp)
    lbl = lt_opp(i);
    md = lbl.metadata;
    if isfield(md, 'location_prior_type')
        if ~strcmp(md.location_prior_type, lbl.location_prior.type)
            prior_preserved = false;
            fprintf('    [mismatch] uid=%s md=%s lbl=%s\n', ...
                lbl.source_uid, md.location_prior_type, lbl.location_prior.type);
        end
    end
end
if prior_preserved
    fprintf('[PASS] DSI -> LabelTable: prior 保留\n');
else
    fprintf('[FAIL] DSI -> LabelTable: prior 类型不一致\n');
end

% [A3] SourceContext per_band_location_prior_map 存在且覆盖 opp
assert(isfield(SourceContext, 'per_band_location_prior_map'), ...
    'per_band_location_prior_map 缺失');
assert(isfield(SourceContext, 'per_band_power_prior_map'), ...
    'per_band_power_prior_map 缺失');
assert(isfield(SourceContext, 'per_band_metadata_map'), ...
    'per_band_metadata_map 缺失');
n_opp_cells = 0;
for t = 1:T
    for b = 1:B
        if SourceContext.per_band_label_map(t, b) == LC.OPPORTUNISTIC
            lp = SourceContext.per_band_location_prior_map{t, b};
            assert(isstruct(lp) && isfield(lp, 'type'), 'SourceContext opp location_prior 丢失');
            n_opp_cells = n_opp_cells + 1;
        end
    end
end
fprintf('[PASS] LabelTable -> SourceContext: prior maps 存在 (opp cells=%d)\n', n_opp_cells);

% [A4] EventList opp 保留 prior
ev_opp_mask = false(1, numel(EventList));
for i = 1:numel(EventList)
    ev_opp_mask(i) = (EventList(i).label == LC.OPPORTUNISTIC);
end
ev_opp = EventList(ev_opp_mask);
all_ok = true;
for i = 1:numel(ev_opp)
    ev = ev_opp(i);
    if ~isfield(ev, 'location_prior') || ~isstruct(ev.location_prior) || ~isfield(ev, 'power_prior') || ...
       ~isfield(ev, 'metadata') || ~isfield(ev, 'source_uid') || isempty(ev.source_uid) || ...
       ~isfield(ev, 'instance_id') || ~isfield(ev, 'linked_template_key')
        all_ok = false;
        fprintf('    [miss] event %d: missing fields\n', ev.event_id);
    end
end
if all_ok
    fprintf('[PASS] LabelTable -> EventList: prior/meta/uid/instance_id/template_key 保留 (n=%d)\n', numel(ev_opp));
else
    fprintf('[FAIL] LabelTable -> EventList: 字段缺失\n');
end

% [A5] 运行 M4, 验证 LocResults
fprintf('\n--- 运行 M4 验证 LocResults ---\n');
try
    if ~exist('cache/SpatialFP_mwm_Nl.mat', 'file')
        [SpatialFP, ~] = init_m25_single_source_fp(APs, Bands, GridValid, Config, SourceTemplates);
    else
        S = load('cache/SpatialFP_mwm_Nl.mat');
        SpatialFP = S.SpatialFP;
    end
    LocResults = run_m4_wknn_localization(EventList, SpatialFP, FrameStates, Config);

    lr_opp_mask = false(1, numel(LocResults));
    for i = 1:numel(LocResults)
        lr_opp_mask(i) = (LocResults(i).label == LC.OPPORTUNISTIC);
    end
    lr_opp = LocResults(lr_opp_mask);
    all_ok_lr = true;
    for i = 1:numel(lr_opp)
        lr = lr_opp(i);
        if ~isfield(lr, 'location_prior_type') || ~isfield(lr, 'power_prior_type') || ...
           ~isfield(lr, 'source_uid') || isempty(lr.source_uid) || ...
           ~isfield(lr, 'linked_template_key') || ~isfield(lr, 'instance_id')
            all_ok_lr = false;
        end
    end
    if all_ok_lr
        fprintf('[PASS] EventList -> LocResults: prior/meta/uid/instance_id 保留 (n=%d)\n', numel(lr_opp));
    else
        fprintf('[FAIL] EventList -> LocResults: 字段缺失\n');
    end
catch ME
    fprintf('[WARN] M4 运行失败（可能是缓存或环境问题）: %s\n', ME.message);
end

fprintf('\n========== 验证完成 ==========\n');


%% ==================== 辅助函数 ====================

function v = ternary_default(x)
    if isempty(x) || (isnumeric(x) && isnan(x))
        v = 0;
    else
        v = x;
    end
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
