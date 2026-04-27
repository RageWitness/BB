%% CB_DEMO_AFTER_MAIN
% Run this script after main_simulation has finished.
% It copies SpatialFP, calibrates only the copied pending library, compares
% localization before/after, and plots per-band opportunistic CB-GP updates.

fprintf('\n============================================\n');
fprintf('  CB AP-profile demo after main_simulation\n');
fprintf('============================================\n\n');

required_vars = {'SpatialFP','SourceContext','Y_dBm_all','APs','Config','EventList','FrameStates'};
for iv = 1:numel(required_vars)
    if ~exist(required_vars{iv}, 'var')
        error('[CB-demo] Missing workspace variable: %s. Run main_simulation first.', required_vars{iv});
    end
end

if ~exist('GridValid', 'var')
    GridValid.xy = SpatialFP.grid_xy;
    GridValid.Nvalid = size(SpatialFP.grid_xy, 1);
end

%% User-facing settings
cb_source_filter = 'opportunistic';     % opportunistic | calibration | all_power_known
cb_train_fraction = 0.80;
cb_min_validation_events = 3;
cb_rng_seed = 20260426;
cb_local_probe_enable = true;
cb_local_probe_radius_m = 60;
cb_local_probe_max_per_band = 30;
cb_local_probe_tx_power_dBm = 20;

% Keep the first run reasonably light. Increase these if you want a stronger
% but slower global AP-profile update.
Config_cb = Config;
Config_cb.calib_ap_profile.enable = true;
Config_cb.calib_ap_profile.power_known_only = true;
Config_cb.calib_ap_profile.ap_pos_xy = APs.pos_xy;
Config_cb.calib_ap_profile.max_offline_train_points = 220;
Config_cb.calib_ap_profile.max_train_points_exact_gp = 330;
Config_cb.calib_ap_profile.prediction_block_size = 1200;
Config_cb.calib_ap_profile.sigma_f_dB = 8;
Config_cb.calib_ap_profile.ell_x_m = 45;
Config_cb.calib_ap_profile.ell_y_m = 45;
Config_cb.calib_ap_profile.sigma0_dB = 8;
Config_cb.calib_ap_profile.sigma_meas_dB = 1.5;
Config_cb.calib_ap_profile.gamma_region_var = 1.0;
Config_cb.calib_ap_profile.loc_rmse_threshold_m = 13;
Config_cb.calib_ap_profile.enable_fast_search = true;
Config_cb.calib_ap_profile.fast_ell_grid_m = [20, 40, 80, 120];
Config_cb.calib_ap_profile.fast_trust_ratio_grid = [2, 4, 8, 12];
Config_cb.calib_ap_profile.fast_refine_enable = true;
Config_cb.calib_ap_profile.enable_ga_feedback = false;
Config_cb.calib_ap_profile.ga_population = 8;
Config_cb.calib_ap_profile.ga_generations = 4;
Config_cb.calib_ap_profile.verbose = true;

rng(cb_rng_seed);
LC = source_label_constants();

%% Copy library and build CB event pool
SpatialFP_original = SpatialFP;
SpatialFP_working = SpatialFP_original;

CBEvents_all = CB_build_calibration_events_from_source_context( ...
    SourceContext, Y_dBm_all, SpatialFP_working, Config_cb);
CBEvents_all = filter_cb_events(CBEvents_all, cb_source_filter, LC);
OppPriorEvents_all = collect_opportunistic_prior_events(SourceContext, LC);

fprintf('[CB-demo] CB event pool after filter=%s: %d\n', ...
    cb_source_filter, numel(CBEvents_all));
fprintf('[CB-demo] all opportunistic prior events for plotting: %d\n', numel(OppPriorEvents_all));
print_cb_event_summary(CBEvents_all, SpatialFP_working.B);

if isempty(CBEvents_all)
    error('[CB-demo] No usable power-known events for CB calibration. Check opportunistic power_prior.exact probability.');
end

[CBEvents_train, ValidationEventList] = split_cb_train_validation( ...
    CBEvents_all, EventList, cb_train_fraction, cb_min_validation_events);

fprintf('[CB-demo] train CB events=%d, validation EventList=%d\n', ...
    numel(CBEvents_train), numel(ValidationEventList));

%% Build AP profiles and calibrate copied library
[APProfile0, SpatialFP_init] = CB_build_ap_profiles_from_offline_library( ...
    SpatialFP_working, SpatialFP_working.grid_xy, Config_cb);

CBSearchResult = struct('status', 'not_run');
Config_cb_final = Config_cb;
cfg_cb_now = CB_calib_ap_profile_defaults(Config_cb);
if ~isempty(ValidationEventList) && cfg_cb_now.enable_fast_search
    [Config_cb_final, CBSearchResult] = CB_search_hyperparams_fast( ...
        SpatialFP_init, CBEvents_train, ValidationEventList, FrameStates, Config_cb);
    [APProfile_pending, SpatialFP_pending, CBDiagnostics] = ...
        CB_update_ap_profiles_with_calibration_sources( ...
            APProfile0, SpatialFP_init, CBEvents_train, Config_cb_final);
elseif isempty(ValidationEventList) || ~cfg_cb_now.enable_ga_feedback
    [APProfile_pending, SpatialFP_pending, CBDiagnostics] = ...
        CB_update_ap_profiles_with_calibration_sources( ...
            APProfile0, SpatialFP_init, CBEvents_train, Config_cb);
else
    [APProfile_pending, SpatialFP_pending, CBDiagnostics] = ...
        CB_run_update_with_loc_feedback( ...
            APProfile0, SpatialFP_init, CBEvents_train, ValidationEventList, FrameStates, Config_cb);
end
CBDiagnostics.fast_search = CBSearchResult;
Config_cb = Config_cb_final;

fprintf('\n--- CB result ---\n');
fprintf('  status: %s\n', CBDiagnostics.status);
if isfield(CBDiagnostics, 'fast_search') && isfield(CBDiagnostics.fast_search, 'status')
    fprintf('  fast search: %s\n', CBDiagnostics.fast_search.status);
    if strcmp(CBDiagnostics.fast_search.status, 'ok')
        fprintf('  best ell=%.2f m, trust_ratio=%.2f, validation RMSE=%.3f m\n', ...
            CBDiagnostics.fast_search.best_ell_m, ...
            CBDiagnostics.fast_search.best_trust_ratio, ...
            CBDiagnostics.fast_search.best_rmse_m);
    end
end
if isfield(CBDiagnostics, 'feedback') && isfield(CBDiagnostics.feedback, 'enabled')
    fprintf('  feedback enabled: %d\n', CBDiagnostics.feedback.enabled);
    if isfield(CBDiagnostics.feedback, 'initial_loc_rmse_m')
        fprintf('  validation RMSE initial: %.3f m\n', CBDiagnostics.feedback.initial_loc_rmse_m);
    end
    if isfield(CBDiagnostics.feedback, 'best_loc_rmse_m')
        fprintf('  validation RMSE best: %.3f m\n', CBDiagnostics.feedback.best_loc_rmse_m);
    end
end

%% Localization before/after using the same EventList
fprintf('\n--- M4 before/after copied library calibration ---\n');
LocResults_before_CB = run_m4_wknn_localization(EventList, SpatialFP_original, FrameStates, Config_cb);
LocResults_after_CB  = run_m4_wknn_localization(EventList, SpatialFP_pending, FrameStates, Config_cb);
print_loc_compare(LocResults_before_CB, LocResults_after_CB);

TargetEventList_CB = filter_eventlist_by_label(EventList, LC.TARGET);
fprintf('\n--- M4 target-only before/after copied library calibration ---\n');
fprintf('  target events: %d / %d\n', numel(TargetEventList_CB), numel(EventList));
LocResults_target_before_CB = run_m4_wknn_localization( ...
    TargetEventList_CB, SpatialFP_original, FrameStates, Config_cb);
LocResults_target_after_CB = run_m4_wknn_localization( ...
    TargetEventList_CB, SpatialFP_pending, FrameStates, Config_cb);
print_loc_compare(LocResults_target_before_CB, LocResults_target_after_CB);

%% Plots
plot_cb_delta_maps(SpatialFP_original, SpatialFP_pending, CBEvents_train, APs, GridValid);
plot_cb_event_map_by_band(CBEvents_train, SpatialFP_original, APs, GridValid);
plot_cb_per_ap_delta(CBDiagnostics, SpatialFP_original);
plot_cb_loc_compare(LocResults_before_CB, LocResults_after_CB, GridValid, APs);
plot_cb_target_loc_with_opp( ...
    LocResults_target_before_CB, LocResults_target_after_CB, CBEvents_train, GridValid, APs);
plot_cb_target_after_by_band(LocResults_target_after_CB, CBEvents_train, SpatialFP_pending, GridValid, APs);

if cb_local_probe_enable
    if exist('Bands', 'var') && exist('ChannelState', 'var')
        [CBLocalProbeEvents, CBLocalProbeFrameStates] = build_cb_local_probe_events( ...
            CBEvents_train, GridValid, APs, Bands, ChannelState, Config_cb, ...
            cb_local_probe_radius_m, cb_local_probe_max_per_band, cb_local_probe_tx_power_dBm);

        fprintf('\n--- CB local probe target before/after around used calibration sources ---\n');
        fprintf('  local probe events: %d, radius=%.1f m\n', ...
            numel(CBLocalProbeEvents), cb_local_probe_radius_m);

        LocResults_local_probe_before_CB = run_m4_wknn_localization( ...
            CBLocalProbeEvents, SpatialFP_original, CBLocalProbeFrameStates, Config_cb);
        LocResults_local_probe_after_CB = run_m4_wknn_localization( ...
            CBLocalProbeEvents, SpatialFP_pending, CBLocalProbeFrameStates, Config_cb);

        print_loc_compare(LocResults_local_probe_before_CB, LocResults_local_probe_after_CB);
        plot_cb_local_probe_compare_by_band( ...
            LocResults_local_probe_before_CB, LocResults_local_probe_after_CB, ...
            CBEvents_train, SpatialFP_pending, GridValid, APs, cb_local_probe_radius_m);
    else
        fprintf('[CB-demo] skip local probe: Bands or ChannelState missing in workspace\n');
    end
end

fprintf('\n[CB-demo] Done. Pending library is in variable SpatialFP_pending.\n');
fprintf('[CB-demo] Original library remains in SpatialFP_original / SpatialFP.\n');


%% ==================== local functions ====================

function events = filter_cb_events(events, mode, LC)
    if isempty(events), return; end
    keep = false(1, numel(events));
    for i = 1:numel(events)
        switch mode
            case 'opportunistic'
                keep(i) = isfield(events(i), 'label') && events(i).label == LC.OPPORTUNISTIC;
            case 'calibration'
                keep(i) = isfield(events(i), 'label') && ...
                    (events(i).label == LC.PERSISTENT_CAL || events(i).label == LC.BROADBAND_CAL);
            case 'all_power_known'
                keep(i) = true;
            otherwise
                error('[CB-demo] Unknown cb_source_filter=%s', mode);
        end
    end
    events = events(keep);
end


function events = collect_opportunistic_prior_events(SourceContext, LC)
    events = struct([]);
    if ~isfield(SourceContext, 'label_table') || isempty(SourceContext.label_table)
        return;
    end
    n = 0;
    for i = 1:numel(SourceContext.label_table)
        lbl = SourceContext.label_table(i);
        if ~isfield(lbl, 'label') || lbl.label ~= LC.OPPORTUNISTIC
            continue;
        end
        if ~isfield(lbl, 'location_prior') || ~isstruct(lbl.location_prior)
            continue;
        end
        ev = struct();
        ev.source_uid = get_label_field(lbl, 'source_uid', sprintf('opp_%d', i));
        ev.band_id = get_label_field(lbl, 'band_id', NaN);
        ev.label = lbl.label;
        ev.location_prior = lbl.location_prior;
        ev.power_prior = get_label_field(lbl, 'power_prior', struct('type', 'none', 'value', []));
        ev.start_frame = get_label_field(lbl, 'start_frame', NaN);
        ev.end_frame = get_label_field(lbl, 'end_frame', NaN);
        ev.metadata = get_label_field(lbl, 'metadata', struct());
        n = n + 1;
        if n == 1
            events = ev;
        else
            events(n) = ev; %#ok<AGROW>
        end
    end
end


function v = get_label_field(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end


function print_cb_event_summary(events, B)
    fprintf('\n--- CB event summary ---\n');
    if isempty(events)
        fprintf('  empty\n');
        return;
    end
    for b = 1:B
        idx = find([events.band_id] == b);
        if isempty(idx)
            fprintf('  band %d: 0\n', b);
            continue;
        end
        types = cell(1, numel(idx));
        for k = 1:numel(idx)
            types{k} = events(idx(k)).location_prior.type;
        end
        fprintf('  band %d: total=%d, exact=%d, region=%d, gaussian=%d, trajectory=%d\n', ...
            b, numel(idx), sum(strcmp(types,'exact')), sum(strcmp(types,'region')), ...
            sum(strcmp(types,'gaussian')), sum(strcmp(types,'trajectory')));
    end
end


function [train_events, val_eventlist] = split_cb_train_validation(cb_events, EventList, train_fraction, min_val)
    train_events = cb_events;
    val_eventlist = [];
    if numel(cb_events) < max(2, min_val + 1) || isempty(EventList)
        return;
    end

    uids = {cb_events.source_uid};
    order = randperm(numel(uids));
    n_train = max(1, floor(train_fraction * numel(uids)));
    n_train = min(n_train, numel(uids) - min_val);
    if n_train < 1
        return;
    end

    train_uid = uids(order(1:n_train));
    val_uid = uids(order(n_train+1:end));
    train_mask = ismember(uids, train_uid);
    train_events = cb_events(train_mask);

    if isempty(val_uid)
        return;
    end

    if isfield(EventList, 'source_uid')
        ev_uids = {EventList.source_uid};
        val_mask = ismember(ev_uids, val_uid);
        val_eventlist = EventList(val_mask);
    end
end


function out = filter_eventlist_by_label(EventList, label_value)
    out = [];
    if isempty(EventList)
        return;
    end
    keep = false(1, numel(EventList));
    for i = 1:numel(EventList)
        keep(i) = isfield(EventList(i), 'label') && EventList(i).label == label_value;
    end
    out = EventList(keep);
end


function print_loc_compare(before, after)
    eb = collect_errors(before);
    ea = collect_errors(after);
    fprintf('  before: N=%d, mean=%.3f, med=%.3f, rmse=%.3f\n', ...
        numel(eb), mean_or_nan(eb), median_or_nan(eb), rmse_or_nan(eb));
    fprintf('  after : N=%d, mean=%.3f, med=%.3f, rmse=%.3f\n', ...
        numel(ea), mean_or_nan(ea), median_or_nan(ea), rmse_or_nan(ea));
    if ~isempty(eb) && ~isempty(ea)
        n = min(numel(eb), numel(ea));
        fprintf('  paired first-%d delta mean(after-before)=%.3f m\n', n, mean(ea(1:n) - eb(1:n)));
    end
end


function err = collect_errors(R)
    if isempty(R) || ~isfield(R, 'loc_error')
        err = [];
        return;
    end
    err = [R.loc_error];
    err = err(isfinite(err));
end


function v = mean_or_nan(x)
    if isempty(x), v = NaN; else, v = mean(x); end
end


function v = median_or_nan(x)
    if isempty(x), v = NaN; else, v = median(x); end
end


function v = rmse_or_nan(x)
    if isempty(x), v = NaN; else, v = sqrt(mean(x.^2)); end
end


function plot_cb_target_loc_with_opp(before, after, opp_events, GridValid, APs)
    figure('Name', 'CB target localization with opportunistic priors', ...
        'Position', [60 70 1500 620]);

    subplot(1,2,1);
    plot_target_loc_panel_with_opp(before, opp_events, GridValid, APs, 'Target localization before CB');

    subplot(1,2,2);
    plot_target_loc_panel_with_opp(after, opp_events, GridValid, APs, 'Target localization after CB pending library');

    sgtitle('All target localization before/after CB, with opportunistic prior sources');
end


function plot_cb_target_after_by_band(LocResults_after, used_events, SpatialFP, GridValid, APs)
    B = SpatialFP.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);

    figure('Name', 'CB target after calibration by band', ...
        'Position', [40 70 1600 780]);

    for b = 1:B
        subplot(n_row, n_col, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', ...
            'Color', [0.93 0.93 0.93], 'MarkerSize', 2);
        hold on;

        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', ...
            'MarkerSize', 9, 'MarkerFaceColor', 'r');

        plot_band_target_results(LocResults_after, b);
        overlay_opportunistic_priors(filter_events_by_band(used_events, b), true);

        axis equal tight;
        grid on;
        title(sprintf('Band %d after CB', b));
        xlabel('X (m)');
        ylabel('Y (m)');

        n_target = count_loc_band(LocResults_after, b);
        n_used = count_event_band(used_events, b);
        text(0.02, 0.98, sprintf('target=%d, used opp=%d', n_target, n_used), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.65]);
    end

    sgtitle('CB 标校后：各单频带 target 定位结果与实际使用的机遇标校源');
end


function plot_band_target_results(LocResults, band_id)
    if isempty(LocResults), return; end
    for k = 1:numel(LocResults)
        lr = LocResults(k);
        if lr.band_id ~= band_id
            continue;
        end
        if isempty(lr.true_pos_xy) || any(~isfinite(lr.true_pos_xy)) || isempty(lr.est_pos_xy)
            continue;
        end
        plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', ...
            'Color', [0.9 0.15 0.15], 'MarkerFaceColor', [0.9 0.15 0.15], ...
            'MarkerSize', 5);
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
            'Color', [0.55 0 0], 'MarkerSize', 7, 'LineWidth', 1.3);
        plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
             [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
            '-', 'Color', [0.85 0.45 0.45], 'LineWidth', 0.8);
        text(lr.true_pos_xy(1) + 2, lr.true_pos_xy(2) + 2, ...
            sprintf('T%d', k), 'Color', [0.8 0 0], 'FontSize', 7);
    end
end


function out = filter_events_by_band(events, band_id)
    out = [];
    if isempty(events), return; end
    keep = [events.band_id] == band_id;
    out = events(keep);
end


function n = count_loc_band(LocResults, band_id)
    if isempty(LocResults)
        n = 0;
    else
        n = sum([LocResults.band_id] == band_id);
    end
end


function n = count_event_band(events, band_id)
    if isempty(events)
        n = 0;
    else
        n = sum([events.band_id] == band_id);
    end
end


function [ProbeEvents, ProbeFrameStates] = build_cb_local_probe_events( ...
    used_events, GridValid, APs, Bands, ChannelState, Config, radius_m, max_per_band, tx_power_dBm)

    LC = source_label_constants();
    B = Config.m0.num_bands;
    M = APs.num;
    ProbeEvents = [];
    probe_pos_by_event = [];
    probe_band_by_event = [];
    eid = 0;

    for b = 1:B
        used_b = filter_events_by_band(used_events, b);
        if isempty(used_b)
            continue;
        end

        mask = false(GridValid.Nvalid, 1);
        for i = 1:numel(used_b)
            xy0 = representative_prior_xy(used_b(i).location_prior);
            if any(~isfinite(xy0))
                continue;
            end
            d = sqrt(sum((GridValid.xy - xy0).^2, 2));
            mask = mask | (d <= radius_m);
        end

        idx = find(mask);
        if isempty(idx)
            continue;
        end
        if numel(idx) > max_per_band
            pick = round(linspace(1, numel(idx), max_per_band));
            idx = idx(pick);
        end

        for ii = 1:numel(idx)
            pos_xy = GridValid.xy(idx(ii), :);
            active_state = make_probe_active_state(b, pos_xy, tx_power_dBm);
            BandObs = generate_band_obs_m1(active_state, APs, Bands, ChannelState, Config, b);

            eid = eid + 1;
            ev = struct();
            ev.event_id = eid;
            ev.band_id = b;
            ev.t_start = eid;
            ev.t_end = eid;
            ev.time_range = [eid, eid];
            ev.duration = 1;
            ev.obs_segment_dBm = BandObs.Y_dBm(:);
            ev.obs_segment_lin = BandObs.Y_lin(:);
            ev.label = LC.TARGET;
            ev.type_hat = 'target';
            ev.route_action = 'localize_only';
            ev.linked_template_key = sprintf('cb_local_probe_b%d_%03d', b, ii);
            ev.source_uid = ev.linked_template_key;
            ev.hold_reason = '';
            ev.upgrade_hint = 'none';
            ev.location_prior = struct('type', 'none', 'value', []);
            ev.power_prior = struct('type', 'exact', 'value', tx_power_dBm);
            ev.metadata = struct('source_type_name', 'target', 'instance_id', eid, ...
                'is_cb_local_probe', true, 'true_pos_xy', pos_xy);
            ev.source_type_name = 'target';
            ev.location_prior_type = 'none';
            ev.power_prior_type = 'exact';
            ev.instance_id = eid;
            ev.is_calibration_source = false;
            ev.is_target_source = true;
            ev.is_opportunistic_source = false;
            ev.power_level_est = mean(BandObs.Y_dBm);
            ev.power_stability_est = 0;
            ev.score_trusted = 0;
            ev.score_prior_pos = 0;
            ev.score_prior_time = 0;
            ev.score_target = 1;
            ev.band_coverage_vec = zeros(1, B);
            ev.band_coverage_vec(b) = 1;
            ev.n_valid_ap = M;

            if isempty(ProbeEvents)
                ProbeEvents = ev;
            else
                ProbeEvents(end+1) = ev; %#ok<AGROW>
            end
            probe_pos_by_event(eid, :) = pos_xy; %#ok<AGROW>
            probe_band_by_event(eid, 1) = b; %#ok<AGROW>
        end
    end

    ProbeFrameStates = build_probe_frame_states(probe_pos_by_event, probe_band_by_event, B, tx_power_dBm);
end


function active_state = make_probe_active_state(b, pos_xy, tx_power_dBm)
    active_state = struct();
    active_state.has_source = true;
    active_state.source_type = 'target';
    active_state.source_subtype = 'target';
    active_state.band_id = b;
    active_state.band_list = b;
    active_state.true_position = pos_xy;
    active_state.true_pos_xy = pos_xy;
    active_state.true_tx_power = tx_power_dBm;
    active_state.tx_power_dBm = tx_power_dBm;
    active_state.tx_power_by_band = [];
    active_state.template_key = sprintf('cb_local_probe_b%d', b);
    active_state.instance_id = 0;
    active_state.location_prior = struct('type', 'none', 'value', []);
    active_state.power_prior = struct('type', 'exact', 'value', tx_power_dBm);
    active_state.timestamp = 1;
    active_state.is_calibrator = false;
    active_state.is_localization_target = true;
end


function FrameStates = build_probe_frame_states(pos_by_event, band_by_event, B, tx_power_dBm)
    T = size(pos_by_event, 1);
    FrameStates = cell(1, T);
    for t = 1:T
        fs = struct();
        fs.frame_id = t;
        fs.time_sec = t;
        for b = 1:B
            fs.active_per_band(b) = empty_probe_event(b, t);
        end
        b0 = band_by_event(t);
        fs.active_per_band(b0).has_source = true;
        fs.active_per_band(b0).source_type = 'target';
        fs.active_per_band(b0).source_subtype = 'target';
        fs.active_per_band(b0).band_id = b0;
        fs.active_per_band(b0).band_list = b0;
        fs.active_per_band(b0).true_position = pos_by_event(t, :);
        fs.active_per_band(b0).true_pos_xy = pos_by_event(t, :);
        fs.active_per_band(b0).true_tx_power = tx_power_dBm;
        fs.active_per_band(b0).tx_power_dBm = tx_power_dBm;
        fs.active_per_band(b0).template_key = sprintf('cb_local_probe_b%d_%03d', b0, t);
        fs.active_per_band(b0).instance_id = t;
        fs.active_per_band(b0).is_localization_target = true;
        fs.active_source_count_unique = 1;
        fs.active_band_count = 1;
        FrameStates{t} = fs;
    end
end


function apb = empty_probe_event(b, t)
    apb = struct();
    apb.has_source = false;
    apb.source_type = '';
    apb.source_subtype = '';
    apb.band_id = b;
    apb.band_list = b;
    apb.true_position = [0 0];
    apb.true_pos_xy = [0 0];
    apb.true_tx_power = 0;
    apb.tx_power_dBm = 0;
    apb.tx_power_by_band = [];
    apb.location_prior = struct('type', 'none', 'value', []);
    apb.power_prior = struct('type', 'none', 'value', []);
    apb.timestamp = t;
    apb.template_key = '';
    apb.instance_id = 0;
    apb.is_calibrator = false;
    apb.is_localization_target = false;
end


function xy = representative_prior_xy(lp)
    switch lp.type
        case 'exact'
            xy = prior_exact_xy(lp);
        case 'region'
            bbox = prior_bbox(lp);
            xy = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
        case 'gaussian'
            [mu, ~] = prior_gaussian(lp);
            xy = mu;
        case 'trajectory'
            path = prior_trajectory_path(lp);
            if isempty(path)
                xy = [NaN NaN];
            else
                xy = path(round(size(path, 1) / 2), :);
            end
        otherwise
            xy = [NaN NaN];
    end
end


function plot_cb_local_probe_compare_by_band(before, after, used_events, SpatialFP, GridValid, APs, radius_m)
    B = SpatialFP.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);

    figure('Name', 'CB local probe before after by band', ...
        'Position', [40 80 1650 800]);

    for b = 1:B
        subplot(n_row, n_col, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', ...
            'Color', [0.94 0.94 0.94], 'MarkerSize', 2);
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', ...
            'MarkerSize', 9, 'MarkerFaceColor', 'r');

        plot_probe_result_pair(before, after, b);
        used_b = filter_events_by_band(used_events, b);
        overlay_probe_radius(used_b, radius_m);
        overlay_opportunistic_priors(used_b, true);

        axis equal tight;
        grid on;
        title(sprintf('Band %d local probes around used opp', b));
        xlabel('X (m)');
        ylabel('Y (m)');
        text(0.02, 0.98, sprintf('probe=%d, used opp=%d', ...
            count_loc_band(after, b), count_event_band(used_events, b)), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.65]);
    end

    sgtitle('60m local synthetic target probes: original library vs CB pending library');
end


function plot_probe_result_pair(before, after, band_id)
    if isempty(after), return; end
    for k = 1:numel(after)
        lr_a = after(k);
        if lr_a.band_id ~= band_id
            continue;
        end
        true_xy = lr_a.true_pos_xy;
        if isempty(true_xy) || any(~isfinite(true_xy))
            continue;
        end
        plot(true_xy(1), true_xy(2), 'o', ...
            'Color', [0.1 0.1 0.1], 'MarkerFaceColor', [0.1 0.1 0.1], ...
            'MarkerSize', 4);

        before_idx = find_matching_loc_result(before, lr_a.event_id);
        if before_idx > 0
            lr_b = before(before_idx);
            plot(lr_b.est_pos_xy(1), lr_b.est_pos_xy(2), 'x', ...
                'Color', [0.75 0.2 0.2], 'MarkerSize', 7, 'LineWidth', 1.2);
            plot([true_xy(1), lr_b.est_pos_xy(1)], [true_xy(2), lr_b.est_pos_xy(2)], ...
                '-', 'Color', [0.9 0.55 0.55], 'LineWidth', 0.7);
        end

        plot(lr_a.est_pos_xy(1), lr_a.est_pos_xy(2), '+', ...
            'Color', [0.05 0.45 0.9], 'MarkerSize', 7, 'LineWidth', 1.3);
        plot([true_xy(1), lr_a.est_pos_xy(1)], [true_xy(2), lr_a.est_pos_xy(2)], ...
            '-', 'Color', [0.45 0.7 0.95], 'LineWidth', 0.7);
    end
end


function idx = find_matching_loc_result(results, event_id)
    idx = 0;
    if isempty(results), return; end
    ids = [results.event_id];
    hit = find(ids == event_id, 1);
    if ~isempty(hit), idx = hit; end
end


function overlay_probe_radius(events, radius_m)
    if isempty(events), return; end
    th = linspace(0, 2*pi, 120);
    for i = 1:numel(events)
        xy = representative_prior_xy(events(i).location_prior);
        if any(~isfinite(xy)), continue; end
        plot(xy(1) + radius_m*cos(th), xy(2) + radius_m*sin(th), ...
            '--', 'Color', [0.25 0.55 0.25], 'LineWidth', 0.8);
    end
end


function plot_target_loc_panel_with_opp(LocResults, opp_events, GridValid, APs, ttl)
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.93 0.93 0.93], 'MarkerSize', 2);
    hold on;

    h_ap = plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', ...
        'MarkerSize', 9, 'MarkerFaceColor', 'r');

    h_t_true = plot(NaN, NaN, 'o', 'Color', [0.9 0.15 0.15], ...
        'MarkerFaceColor', [0.9 0.15 0.15], 'MarkerSize', 6);
    h_t_est = plot(NaN, NaN, 'x', 'Color', [0.55 0 0], ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    h_line = plot(NaN, NaN, '-', 'Color', [0.85 0.45 0.45], 'LineWidth', 0.8);

    if ~isempty(LocResults)
        for k = 1:numel(LocResults)
            lr = LocResults(k);
            if isempty(lr.true_pos_xy) || any(~isfinite(lr.true_pos_xy)) || isempty(lr.est_pos_xy)
                continue;
            end
            plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', ...
                'Color', [0.9 0.15 0.15], 'MarkerFaceColor', [0.9 0.15 0.15], 'MarkerSize', 5);
            plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
                'Color', [0.55 0 0], 'MarkerSize', 7, 'LineWidth', 1.3);
            plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
                 [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
                '-', 'Color', [0.85 0.45 0.45], 'LineWidth', 0.8);
            text(lr.true_pos_xy(1) + 2, lr.true_pos_xy(2) + 2, ...
                sprintf('B%d target', lr.band_id), ...
                'Color', [0.8 0 0], 'FontSize', 7);
        end
    end

    h_opp = overlay_opportunistic_priors(opp_events, true);

    axis equal tight;
    grid on;
    title(ttl);
    xlabel('X (m)');
    ylabel('Y (m)');

    leg_h = gobjects(0);
    leg_s = {};
    leg_h(end+1) = h_ap; %#ok<AGROW>
    leg_s{end+1} = 'AP'; %#ok<AGROW>
    leg_h(end+1) = h_t_true; %#ok<AGROW>
    leg_s{end+1} = 'target true'; %#ok<AGROW>
    leg_h(end+1) = h_t_est; %#ok<AGROW>
    leg_s{end+1} = 'target est'; %#ok<AGROW>
    leg_h(end+1) = h_line; %#ok<AGROW>
    leg_s{end+1} = 'error line'; %#ok<AGROW>
    if ~isempty(h_opp.exact)
        leg_h(end+1) = h_opp.exact; %#ok<AGROW>
        leg_s{end+1} = 'opp exact'; %#ok<AGROW>
    end
    if ~isempty(h_opp.region)
        leg_h(end+1) = h_opp.region; %#ok<AGROW>
        leg_s{end+1} = 'opp region'; %#ok<AGROW>
    end
    if ~isempty(h_opp.gaussian)
        leg_h(end+1) = h_opp.gaussian; %#ok<AGROW>
        leg_s{end+1} = 'opp gaussian'; %#ok<AGROW>
    end
    if ~isempty(h_opp.trajectory)
        leg_h(end+1) = h_opp.trajectory; %#ok<AGROW>
        leg_s{end+1} = 'opp trajectory'; %#ok<AGROW>
    end
    add_safe_legend(leg_h, leg_s);
end


function h = overlay_opportunistic_priors(events, show_text)
    h = struct('exact', [], 'region', [], 'gaussian', [], 'trajectory', []);
    if isempty(events), return; end

    for i = 1:numel(events)
        if ~isfield(events(i), 'location_prior') || ~isstruct(events(i).location_prior)
            continue;
        end
        lp = events(i).location_prior;
        b = events(i).band_id;
        switch lp.type
            case 'exact'
                xy = prior_exact_xy(lp);
                hh = plot(xy(1), xy(2), 'd', 'Color', [0.75 0.45 0], ...
                    'MarkerFaceColor', [1 0.75 0], 'MarkerSize', 6, 'LineWidth', 1.0);
                if isempty(h.exact), h.exact = hh; end
                label_xy = xy;
                label_txt = sprintf('B%d opp exact', b);

            case 'region'
                bbox = prior_bbox(lp);
                rectangle('Position', [bbox(1), bbox(3), bbox(2)-bbox(1), bbox(4)-bbox(3)], ...
                    'EdgeColor', [0 0.55 0], 'LineWidth', 1.2, 'LineStyle', '-');
                if isempty(h.region)
                    h.region = plot(NaN, NaN, '-', 'Color', [0 0.55 0], ...
                        'LineWidth', 1.2, 'LineStyle', '-');
                end
                label_xy = [bbox(1), bbox(4)];
                label_txt = sprintf('B%d opp region', b);

            case 'gaussian'
                [mu, Sigma] = prior_gaussian(lp);
                hh = plot_gaussian_ellipse(mu, Sigma, [0.1 0.2 0.9]);
                plot(mu(1), mu(2), 'x', 'Color', [0.1 0.2 0.9], 'LineWidth', 1.3);
                if isempty(h.gaussian), h.gaussian = hh; end
                label_xy = mu;
                label_txt = sprintf('B%d opp gauss', b);

            case 'trajectory'
                path = prior_trajectory_path(lp);
                if isempty(path), continue; end
                hh = plot(path(:,1), path(:,2), '-', 'Color', [0.55 0 0.75], 'LineWidth', 1.4);
                if isempty(h.trajectory), h.trajectory = hh; end
                label_xy = path(round(size(path,1)/2), :);
                label_txt = sprintf('B%d opp traj', b);

            otherwise
                continue;
        end

        if show_text
            text(label_xy(1) + 2, label_xy(2) + 2, label_txt, ...
                'Color', [0.15 0.15 0.15], 'FontSize', 7);
        end
    end
end


function add_safe_legend(handles_in, labels_in)
    keep = false(1, numel(handles_in));
    for i = 1:numel(handles_in)
        keep(i) = ~isempty(handles_in(i)) && isgraphics(handles_in(i));
    end
    if any(keep)
        legend(handles_in(keep), labels_in(keep), 'Location', 'best');
    end
end


function path = prior_trajectory_path(lp)
    if isfield(lp, 'path_points')
        path = lp.path_points;
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'path_points')
        path = lp.value.path_points;
    else
        path = [];
    end
end


function plot_cb_delta_maps(SpatialFP_old, SpatialFP_new, events, APs, GridValid)
    B = SpatialFP_old.B;
    vals = cell(1, B);
    clim_max = 0;
    for b = 1:B
        D = SpatialFP_new.band(b).F_dBm - SpatialFP_old.band(b).F_dBm;
        vals{b} = sqrt(mean(D.^2, 1));
        clim_max = max(clim_max, prctile(abs(vals{b}(:)), 98));
    end
    if clim_max <= 0 || ~isfinite(clim_max), clim_max = 1; end

    figure('Name', 'CB per-band GP update maps', 'Position', [40 80 1500 420]);
    for b = 1:B
        subplot(1, B, b);
        scatter(GridValid.xy(:,1), GridValid.xy(:,2), 12, vals{b}, 'filled');
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
        overlay_events_for_band(events, b);
        axis equal tight;
        grid on;
        colorbar;
        caxis([0, clim_max]);
        title(sprintf('Band %d RMS update dB', b));
        xlabel('X (m)');
        ylabel('Y (m)');
    end
    sgtitle('CB GP update caused by training opportunistic events');
end


function plot_cb_event_map_by_band(events, SpatialFP, APs, GridValid)
    B = SpatialFP.B;
    figure('Name', 'CB training opportunistic priors by band', 'Position', [80 120 1500 420]);
    for b = 1:B
        subplot(1, B, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 2);
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
        overlay_events_for_band(events, b);
        axis equal tight;
        grid on;
        title(sprintf('Band %d CB events', b));
        xlabel('X (m)');
        ylabel('Y (m)');
    end
    sgtitle('CB training events and location priors');
end


function overlay_events_for_band(events, b)
    if isempty(events), return; end
    for i = 1:numel(events)
        if events(i).band_id ~= b, continue; end
        lp = events(i).location_prior;
        switch lp.type
            case 'exact'
                xy = prior_exact_xy(lp);
                plot(xy(1), xy(2), 'ko', 'MarkerFaceColor', [1 0.75 0], 'MarkerSize', 6);
            case 'region'
                bbox = prior_bbox(lp);
                rectangle('Position', [bbox(1), bbox(3), bbox(2)-bbox(1), bbox(4)-bbox(3)], ...
                    'EdgeColor', [0 0.6 0], 'LineWidth', 1.2, 'LineStyle', '-');
            case 'gaussian'
                [mu, Sigma] = prior_gaussian(lp);
                plot_gaussian_ellipse(mu, Sigma, [0.1 0.2 0.9]);
                plot(mu(1), mu(2), 'x', 'Color', [0.1 0.2 0.9], 'LineWidth', 1.5);
            case 'trajectory'
                if isfield(lp, 'path_points')
                    path = lp.path_points;
                elseif isfield(lp, 'value') && isfield(lp.value, 'path_points')
                    path = lp.value.path_points;
                else
                    path = [];
                end
                if ~isempty(path)
                    plot(path(:,1), path(:,2), '-', 'Color', [0.6 0 0.8], 'LineWidth', 1.5);
                end
        end
    end
end


function xy = prior_exact_xy(lp)
    if isfield(lp, 'xy'), xy = lp.xy(:)';
    else, xy = lp.value(:)'; end
end


function bbox = prior_bbox(lp)
    if isfield(lp, 'bbox')
        bbox = lp.bbox(:)';
    elseif isnumeric(lp.value)
        bbox = lp.value(:)';
    elseif isfield(lp.value, 'bbox')
        bbox = lp.value.bbox(:)';
    else
        bbox = lp.value.region_geometry(:)';
    end
    bbox = [min(bbox(1), bbox(2)), max(bbox(1), bbox(2)), ...
            min(bbox(3), bbox(4)), max(bbox(3), bbox(4))];
end


function [mu, Sigma] = prior_gaussian(lp)
    if isfield(lp, 'mu'), mu = lp.mu(:)';
    else, mu = lp.value.mu(:)'; end
    if isfield(lp, 'Sigma')
        Sigma = lp.Sigma;
    elseif isfield(lp, 'sigma')
        Sigma = eye(2) * lp.sigma^2;
    elseif isfield(lp.value, 'Sigma')
        Sigma = lp.value.Sigma;
    else
        Sigma = eye(2) * lp.value.sigma^2;
    end
end


function h = plot_gaussian_ellipse(mu, Sigma, color)
    th = linspace(0, 2*pi, 80);
    [V, D] = eig(0.5 * (Sigma + Sigma'));
    r = sqrt(max(diag(D), 0));
    pts = [cos(th(:)), sin(th(:))] * diag(r) * V' + mu;
    h = plot(pts(:,1), pts(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end


function plot_cb_per_ap_delta(diag, SpatialFP)
    B = SpatialFP.B;
    M = SpatialFP.M;
    A = nan(M, B);
    for b = 1:B
        for m = 1:M
            if isfield(diag.band(b).ap(m), 'mean_abs_delta')
                A(m, b) = diag.band(b).ap(m).mean_abs_delta;
            end
        end
    end
    figure('Name', 'CB per-AP mean update', 'Position', [100 120 800 420]);
    bar(A);
    grid on;
    xlabel('AP index');
    ylabel('Mean abs update (dB)');
    legend(make_band_labels(B), 'Location', 'best');
    title('CB mean absolute AP-profile update per AP');
end


function labels = make_band_labels(B)
    labels = cell(1, B);
    for b = 1:B
        labels{b} = sprintf('Band %d', b);
    end
end


function plot_cb_loc_compare(before, after, GridValid, APs)
    figure('Name', 'CB localization before after', 'Position', [80 80 1450 560]);
    subplot(1,2,1);
    plot_loc_panel(before, GridValid, APs, 'Before CB');
    subplot(1,2,2);
    plot_loc_panel(after, GridValid, APs, 'After CB pending library');
    sgtitle('M4 localization with original library vs CB pending library');
end


function plot_loc_panel(LocResults, GridValid, APs, ttl)
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.92 0.92 0.92], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
    if ~isempty(LocResults)
        for k = 1:numel(LocResults)
            lr = LocResults(k);
            if isempty(lr.true_pos_xy) || any(~isfinite(lr.true_pos_xy)) || isempty(lr.est_pos_xy)
                continue;
            end
            c = type_color(lr);
            plot(lr.true_pos_xy(1), lr.true_pos_xy(2), 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5);
            plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', 'Color', c * 0.65, 'MarkerSize', 7, 'LineWidth', 1.3);
            line_color = 0.65 * c + 0.35 * [1 1 1];
            plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
                '-', 'Color', line_color, 'LineWidth', 0.8);
        end
    end
    axis equal tight;
    grid on;
    title(ttl);
    xlabel('X (m)');
    ylabel('Y (m)');
end


function c = type_color(lr)
    if isfield(lr, 'is_opportunistic_source') && lr.is_opportunistic_source
        c = [1 0.65 0];
    elseif isfield(lr, 'is_target_source') && lr.is_target_source
        c = [0.9 0.15 0.15];
    else
        c = [0.2 0.2 0.2];
    end
end
