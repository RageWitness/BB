%% CD_DEMO_AFTER_MAIN
% Run this script after main_simulation has finished.
% It copies SpatialFP, applies CD uncertain opportunistic calibration to the
% copied pending library, then compares localization before and after.

fprintf('\n============================================\n');
fprintf('  CD uncertain opportunistic demo after main_simulation\n');
fprintf('============================================\n\n');

required_vars = {'SpatialFP','SourceContext','Y_dBm_all','APs','Config','EventList','FrameStates'};
for iv = 1:numel(required_vars)
    if ~exist(required_vars{iv}, 'var')
        error('[CD-demo] Missing workspace variable: %s. Run main_simulation first.', required_vars{iv});
    end
end

if ~exist('GridValid', 'var')
    GridValid.xy = SpatialFP.grid_xy;
    GridValid.Nvalid = size(SpatialFP.grid_xy, 1);
end

%% User settings
cd_rng_seed = 20260428;
cd_local_probe_enable = true;
cd_local_probe_radius_m = 80;
cd_local_probe_max_per_band = 30;
cd_local_probe_tx_power_dBm = 20;
cd_truth_fp_path = fullfile(pwd, 'cache', 'SpatialFP_lognormal_new.mat');

Config_cd = Config;
Config_cd.cd.enable = true;
Config_cd.cd.source_filter = 'opportunistic';
Config_cd.cd.supported_prior_types = {'region', 'gaussian'};
Config_cd.cd.patch_radius_m = 80;
Config_cd.cd.patch_taper_width_m = 20;
Config_cd.cd.max_base_points = 8;
Config_cd.cd.max_ray_points = 160;
Config_cd.cd.ray_width_m = 9;
Config_cd.cd.blend_old_weight = 1.0;
Config_cd.cd.calibration_weight = 2.0;
Config_cd.cd.max_delta_dB = 25;
Config_cd.cd.verbose = true;

rng(cd_rng_seed);
LC = source_label_constants();

%% Copy library and build CD event pool
SpatialFP_original = SpatialFP;
SpatialFP_working = SpatialFP_original;

AllPowerKnownEvents = CB_build_calibration_events_from_source_context( ...
    SourceContext, Y_dBm_all, SpatialFP_working, Config_cd);
CDEvents_train = filter_cd_events(AllPowerKnownEvents, LC);

fprintf('[CD-demo] power-known event pool: %d\n', numel(AllPowerKnownEvents));
fprintf('[CD-demo] CD uncertain opportunistic events: %d\n', numel(CDEvents_train));
print_cd_event_summary(CDEvents_train, SpatialFP_working.B);

if isempty(CDEvents_train)
    error('[CD-demo] No usable CD events. Need power-known opportunistic sources with region/gaussian location_prior.');
end

%% Calibrate copied library
[SpatialFP_CD, CDResult] = CD_calibrate_spatialfp_opportunistic_uncertain( ...
    SpatialFP_working, CDEvents_train, APs, Config_cd);

fprintf('\n--- CD result ---\n');
fprintf('  status: %s\n', CDResult.status);
fprintf('  valid events: %d / %d\n', CDResult.n_events_valid, CDResult.n_events_total);
for b = 1:SpatialFP_CD.B
    fprintf('  band %d: %s, updated=%d, mean_abs_delta=%.3f dB, max_abs_delta=%.3f dB, n0_median=%.3f\n', ...
        b, CDResult.band(b).status, CDResult.band(b).n_entries_updated, ...
        CDResult.band(b).mean_abs_delta_dB, CDResult.band(b).max_abs_delta_dB, ...
        CDResult.band(b).n0_hat_median);
end

%% Localization before/after using the same EventList
fprintf('\n--- M4 before/after CD copied library calibration ---\n');
LocResults_before_CD = run_m4_wknn_localization(EventList, SpatialFP_original, FrameStates, Config_cd);
LocResults_after_CD  = run_m4_wknn_localization(EventList, SpatialFP_CD, FrameStates, Config_cd);
print_loc_compare(LocResults_before_CD, LocResults_after_CD);

TargetEventList_CD = filter_eventlist_by_label(EventList, LC.TARGET);
fprintf('\n--- M4 target-only before/after CD copied library calibration ---\n');
fprintf('  target events: %d / %d\n', numel(TargetEventList_CD), numel(EventList));
LocResults_target_before_CD = run_m4_wknn_localization( ...
    TargetEventList_CD, SpatialFP_original, FrameStates, Config_cd);
LocResults_target_after_CD = run_m4_wknn_localization( ...
    TargetEventList_CD, SpatialFP_CD, FrameStates, Config_cd);
print_loc_compare(LocResults_target_before_CD, LocResults_target_after_CD);

%% Plots
plot_cd_delta_maps(SpatialFP_original, SpatialFP_CD, CDEvents_train, APs, GridValid);
SpatialFP_truth_CD = load_cd_truth_spatialfp(cd_truth_fp_path);
if ~isempty(SpatialFP_truth_CD)
    plot_cd_truth_error_in_used_regions( ...
        SpatialFP_original, SpatialFP_CD, SpatialFP_truth_CD, ...
        CDEvents_train, APs, GridValid, cd_local_probe_radius_m);
else
    fprintf('[CD-demo] skip truth error map: cannot load %s\n', cd_truth_fp_path);
end
plot_cd_event_map_by_band(CDEvents_train, SpatialFP_original, APs, GridValid);
plot_cd_target_loc_compare( ...
    LocResults_target_before_CD, LocResults_target_after_CD, CDEvents_train, GridValid, APs);
plot_cd_target_after_by_band(LocResults_target_after_CD, CDEvents_train, SpatialFP_CD, GridValid, APs);

if cd_local_probe_enable
    if exist('Bands', 'var') && exist('ChannelState', 'var')
        [CDLocalProbeEvents, CDLocalProbeFrameStates] = build_cd_local_probe_events( ...
            CDEvents_train, GridValid, APs, Bands, ChannelState, Config_cd, ...
            cd_local_probe_radius_m, cd_local_probe_max_per_band, cd_local_probe_tx_power_dBm);

        fprintf('\n--- CD local probe target before/after around used uncertain opportunistic sources ---\n');
        fprintf('  local probe events: %d, radius=%.1f m\n', ...
            numel(CDLocalProbeEvents), cd_local_probe_radius_m);

        LocResults_local_probe_before_CD = run_m4_wknn_localization( ...
            CDLocalProbeEvents, SpatialFP_original, CDLocalProbeFrameStates, Config_cd);
        LocResults_local_probe_after_CD = run_m4_wknn_localization( ...
            CDLocalProbeEvents, SpatialFP_CD, CDLocalProbeFrameStates, Config_cd);

        print_loc_compare(LocResults_local_probe_before_CD, LocResults_local_probe_after_CD);
        print_loc_compare_by_band(LocResults_local_probe_before_CD, LocResults_local_probe_after_CD, SpatialFP_CD.B);
        plot_cd_local_probe_compare_by_band( ...
            LocResults_local_probe_before_CD, LocResults_local_probe_after_CD, ...
            CDEvents_train, SpatialFP_CD, GridValid, APs, cd_local_probe_radius_m);
    else
        fprintf('[CD-demo] skip local probe: Bands or ChannelState missing in workspace\n');
    end
end

fprintf('\n[CD-demo] Done. Pending library is in variable SpatialFP_CD.\n');
fprintf('[CD-demo] Original library remains in SpatialFP_original / SpatialFP.\n');


%% ==================== local functions ====================

function events = filter_cd_events(events, LC)
    if isempty(events), return; end
    keep = false(1, numel(events));
    for i = 1:numel(events)
        if ~isfield(events(i), 'label') || events(i).label ~= LC.OPPORTUNISTIC
            continue;
        end
        if ~isfield(events(i), 'location_prior') || ~isstruct(events(i).location_prior) || ...
                ~isfield(events(i).location_prior, 'type')
            continue;
        end
        keep(i) = any(strcmpi(events(i).location_prior.type, {'region', 'gaussian'}));
    end
    events = events(keep);
end


function print_cd_event_summary(events, B)
    fprintf('\n--- CD event summary ---\n');
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
        fprintf('  band %d: total=%d, region=%d, gaussian=%d\n', ...
            b, numel(idx), sum(strcmpi(types, 'region')), sum(strcmpi(types, 'gaussian')));
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


function print_loc_compare_by_band(before, after, B)
    fprintf('\n--- Local probe before/after by band ---\n');
    for b = 1:B
        eb = collect_errors_by_band(before, b);
        ea = collect_errors_by_band(after, b);
        fprintf('  band %d: N=%d, before mean/rmse=%.3f/%.3f, after mean/rmse=%.3f/%.3f, delta mean=%.3f m\n', ...
            b, numel(ea), mean_or_nan(eb), rmse_or_nan(eb), ...
            mean_or_nan(ea), rmse_or_nan(ea), mean_or_nan(ea) - mean_or_nan(eb));
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


function err = collect_errors_by_band(R, band_id)
    if isempty(R) || ~isfield(R, 'loc_error')
        err = [];
        return;
    end
    keep = [R.band_id] == band_id;
    rr = R(keep);
    if isempty(rr)
        err = [];
    else
        err = [rr.loc_error];
        err = err(isfinite(err));
    end
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


function SpatialFP_truth = load_cd_truth_spatialfp(mat_path)
    SpatialFP_truth = [];
    if ~exist(mat_path, 'file')
        return;
    end

    S = load(mat_path);
    if isfield(S, 'SpatialFP')
        SpatialFP_truth = S.SpatialFP;
        return;
    end

    names = fieldnames(S);
    for i = 1:numel(names)
        v = S.(names{i});
        if isstruct(v) && isfield(v, 'band')
            SpatialFP_truth = v;
            return;
        end
    end
end


function plot_cd_delta_maps(SpatialFP_old, SpatialFP_new, events, APs, GridValid)
    B = SpatialFP_old.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);
    vals = [];
    for b = 1:B
        d = get_band_dbm(SpatialFP_new.band(b)) - get_band_dbm(SpatialFP_old.band(b));
        vals = [vals; abs(d(:))]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
    if isempty(vals)
        clim = [0 1];
    else
        hi = max(vals);
        if hi <= 0, hi = 1; end
        clim = [0 hi];
    end

    figure('Name', 'CD fingerprint delta maps', 'Position', [50 90 1550 520]);
    for b = 1:B
        subplot(n_row, n_col, b);
        Fd = get_band_dbm(SpatialFP_new.band(b)) - get_band_dbm(SpatialFP_old.band(b));
        rms_delta = sqrt(mean(Fd.^2, 1, 'omitnan'));
        scatter(GridValid.xy(:,1), GridValid.xy(:,2), 14, rms_delta(:), 'filled');
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        overlay_opportunistic_priors(filter_events_by_band(events, b), true);
        axis equal tight;
        grid on;
        caxis(clim);
        colorbar;
        title(sprintf('Band %d CD RMS update dB', b));
        xlabel('X (m)');
        ylabel('Y (m)');
    end
    sgtitle('CD calibration delta maps');
end


function plot_cd_truth_error_in_used_regions(SpatialFP_old, SpatialFP_new, SpatialFP_truth, events, APs, GridValid, radius_m)
    B = min([SpatialFP_old.B, SpatialFP_new.B, SpatialFP_truth.B]);
    n_col = B;
    n_row = 2;
    vals = [];
    stat = struct([]);

    for b = 1:B
        F_old = align_band_matrix(get_band_dbm(SpatialFP_old.band(b)), GridValid);
        F_new = align_band_matrix(get_band_dbm(SpatialFP_new.band(b)), GridValid);
        F_true = align_band_matrix(get_band_dbm(SpatialFP_truth.band(b)), GridValid);
        mask = build_cd_used_region_mask(GridValid.xy, filter_events_by_band(events, b), radius_m);

        err_old = rms_ap_error(F_old, F_true);
        err_new = rms_ap_error(F_new, F_true);
        err_old(~mask) = NaN;
        err_new(~mask) = NaN;

        vals = [vals; err_old(isfinite(err_old)); err_new(isfinite(err_new))]; %#ok<AGROW>
        stat(b).mask = mask; %#ok<AGROW>
        stat(b).err_old = err_old;
        stat(b).err_new = err_new;
        stat(b).mean_old = mean_or_nan(err_old(isfinite(err_old)));
        stat(b).mean_new = mean_or_nan(err_new(isfinite(err_new)));
        stat(b).delta_mean = stat(b).mean_new - stat(b).mean_old;
        stat(b).improve_count = sum(isfinite(err_old) & isfinite(err_new) & err_new < err_old);
        stat(b).n_mask = nnz(mask);
    end

    if isempty(vals)
        clim = [0 1];
    else
        hi = percentile_local(vals, 95);
        if ~isfinite(hi) || hi <= 0, hi = max(vals); end
        if ~isfinite(hi) || hi <= 0, hi = 1; end
        clim = [0 hi];
    end

    figure('Name', 'CD truth-library error in used regions', 'Position', [30 70 1700 760]);
    for b = 1:B
        subplot(n_row, n_col, b);
        scatter(GridValid.xy(:,1), GridValid.xy(:,2), 14, stat(b).err_old(:), 'filled');
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        overlay_probe_radius(filter_events_by_band(events, b), radius_m);
        overlay_opportunistic_priors(filter_events_by_band(events, b), true);
        axis equal tight;
        grid on;
        caxis(clim);
        colorbar;
        title(sprintf('B%d before truth err, mean=%.2f dB', b, stat(b).mean_old));
        xlabel('X (m)');
        ylabel('Y (m)');

        subplot(n_row, n_col, B + b);
        scatter(GridValid.xy(:,1), GridValid.xy(:,2), 14, stat(b).err_new(:), 'filled');
        hold on;
        imp = isfinite(stat(b).err_old) & isfinite(stat(b).err_new) & stat(b).err_new < stat(b).err_old;
        plot(GridValid.xy(imp,1), GridValid.xy(imp,2), '.', 'Color', [0.0 0.75 0.2], 'MarkerSize', 6);
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        overlay_probe_radius(filter_events_by_band(events, b), radius_m);
        overlay_opportunistic_priors(filter_events_by_band(events, b), true);
        axis equal tight;
        grid on;
        caxis(clim);
        colorbar;
        title(sprintf('B%d after truth err, mean=%.2f dB, d=%.2f', ...
            b, stat(b).mean_new, stat(b).delta_mean));
        xlabel('X (m)');
        ylabel('Y (m)');
        text(0.02, 0.98, sprintf('inside=%d, improved=%d', stat(b).n_mask, stat(b).improve_count), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.65]);
    end
    sgtitle(sprintf('RF raw dBm RMS error against SpatialFP\\_lognormal\\_new inside used CD regions (radius %.0fm)', radius_m));
end


function mask = build_cd_used_region_mask(grid_xy, events, radius_m)
    mask = false(size(grid_xy, 1), 1);
    if isempty(events), return; end
    for i = 1:numel(events)
        lp = events(i).location_prior;
        local_mask = false(size(mask));
        switch lower(lp.type)
            case 'region'
                bbox = prior_bbox(lp);
                local_mask = grid_xy(:,1) >= bbox(1) & grid_xy(:,1) <= bbox(2) & ...
                    grid_xy(:,2) >= bbox(3) & grid_xy(:,2) <= bbox(4);
                center_xy = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
            case 'gaussian'
                [mu, Sigma] = prior_gaussian(lp);
                S = 0.5 * (Sigma + Sigma') + 1e-9 * eye(2);
                X = grid_xy - mu;
                local_mask = sum((X / S) .* X, 2) <= 5.991;
                center_xy = mu;
            otherwise
                center_xy = representative_prior_xy(lp);
        end
        if all(isfinite(center_xy))
            d = sqrt(sum((grid_xy - center_xy).^2, 2));
            local_mask = local_mask | (d <= radius_m);
        end
        mask = mask | local_mask;
    end
end


function e = rms_ap_error(F, F_ref)
    M = min(size(F, 1), size(F_ref, 1));
    G = min(size(F, 2), size(F_ref, 2));
    D = F(1:M, 1:G) - F_ref(1:M, 1:G);
    e = NaN(G, 1);
    for g = 1:G
        v = D(:, g);
        v = v(isfinite(v));
        if ~isempty(v)
            e(g) = sqrt(mean(v.^2));
        end
    end
end


function F = align_band_matrix(F, GridValid)
    G = size(GridValid.xy, 1);
    if size(F, 2) ~= G && size(F, 1) == G
        F = F';
    end
    if size(F, 2) > G
        F = F(:, 1:G);
        return;
    end
    if size(F, 2) < G
        pad = NaN(size(F, 1), G - size(F, 2));
        F = [F, pad];
    end
end


function q = percentile_local(x, pct)
    x = sort(x(isfinite(x)));
    if isempty(x)
        q = NaN;
        return;
    end
    pos = 1 + (numel(x) - 1) * pct / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (x(hi) - x(lo)) * (pos - lo);
    end
end


function plot_cd_event_map_by_band(events, SpatialFP, APs, GridValid)
    B = SpatialFP.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);
    figure('Name', 'CD used opportunistic priors by band', 'Position', [80 80 1550 650]);
    for b = 1:B
        subplot(n_row, n_col, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.94 0.94 0.94], 'MarkerSize', 2);
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        overlay_opportunistic_priors(filter_events_by_band(events, b), true);
        axis equal tight;
        grid on;
        title(sprintf('Band %d used CD priors', b));
        xlabel('X (m)');
        ylabel('Y (m)');
    end
    sgtitle('CD used uncertain opportunistic sources');
end


function plot_cd_target_loc_compare(before, after, events, GridValid, APs)
    figure('Name', 'CD target localization before after', 'Position', [60 80 1500 620]);
    subplot(1,2,1);
    plot_target_panel(before, events, GridValid, APs, 'Target localization before CD');
    subplot(1,2,2);
    plot_target_panel(after, events, GridValid, APs, 'Target localization after CD pending library');
    sgtitle('Target localization before/after CD');
end


function plot_cd_target_after_by_band(LocResults_after, used_events, SpatialFP, GridValid, APs)
    B = SpatialFP.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);
    figure('Name', 'CD target after calibration by band', 'Position', [40 70 1600 780]);
    for b = 1:B
        subplot(n_row, n_col, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.94 0.94 0.94], 'MarkerSize', 2);
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        plot_band_loc_results(LocResults_after, b);
        overlay_opportunistic_priors(filter_events_by_band(used_events, b), true);
        axis equal tight;
        grid on;
        title(sprintf('Band %d after CD', b));
        xlabel('X (m)');
        ylabel('Y (m)');
        add_cd_target_legend();
        text(0.02, 0.98, sprintf('target=%d, used opp=%d', ...
            count_loc_band(LocResults_after, b), count_event_band(used_events, b)), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.65]);
    end
    sgtitle('CD calibration result by band');
end


function [ProbeEvents, ProbeFrameStates] = build_cd_local_probe_events( ...
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
            ev.linked_template_key = sprintf('cd_local_probe_b%d_%03d', b, ii);
            ev.source_uid = ev.linked_template_key;
            ev.hold_reason = '';
            ev.upgrade_hint = 'none';
            ev.location_prior = struct('type', 'none', 'value', []);
            ev.power_prior = struct('type', 'exact', 'value', tx_power_dBm);
            ev.metadata = struct('source_type_name', 'target', 'instance_id', eid, ...
                'is_cd_local_probe', true, 'true_pos_xy', pos_xy);
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
    active_state.template_key = sprintf('cd_local_probe_b%d', b);
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
        fs.active_per_band(b0).template_key = sprintf('cd_local_probe_b%d_%03d', b0, t);
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


function plot_cd_local_probe_compare_by_band(before, after, used_events, SpatialFP, GridValid, APs, radius_m)
    plot_cd_local_probe_one_figure(before, after, used_events, SpatialFP, GridValid, APs, radius_m, 'before');
    plot_cd_local_probe_one_figure(after, before, used_events, SpatialFP, GridValid, APs, radius_m, 'after');
end


function plot_cd_local_probe_one_figure(primary, other, used_events, SpatialFP, GridValid, APs, radius_m, mode_name)
    B = SpatialFP.B;
    n_col = min(B, 4);
    n_row = ceil(B / n_col);
    if strcmpi(mode_name, 'after')
        fig_name = 'CD local probe after by band';
        main_title = sprintf('%.0fm local synthetic target probes: CD pending library', radius_m);
    else
        fig_name = 'CD local probe before by band';
        main_title = sprintf('%.0fm local synthetic target probes: original library', radius_m);
    end
    figure('Name', fig_name, 'Position', [40 80 1650 800]);
    for b = 1:B
        subplot(n_row, n_col, b);
        plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.94 0.94 0.94], 'MarkerSize', 2);
        hold on;
        plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        [n_probe, n_improved] = plot_probe_result_single(primary, other, b, mode_name);
        used_b = filter_events_by_band(used_events, b);
        overlay_probe_radius(used_b, radius_m);
        overlay_opportunistic_priors(used_b, true);
        axis equal tight;
        grid on;
        title(sprintf('Band %d local probes %s', b, mode_name));
        xlabel('X (m)');
        ylabel('Y (m)');
        add_cd_probe_legend(mode_name);
        text(0.02, 0.98, sprintf('probe=%d, improved=%d, used opp=%d', ...
            n_probe, n_improved, count_event_band(used_events, b)), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.65]);
    end
    sgtitle(main_title);
end


function [n_probe, n_improved] = plot_probe_result_single(primary, other, band_id, mode_name)
    n_probe = 0;
    n_improved = 0;
    if isempty(primary), return; end
    for k = 1:numel(primary)
        lr = primary(k);
        if lr.band_id ~= band_id
            continue;
        end
        true_xy = lr.true_pos_xy;
        if isempty(true_xy) || any(~isfinite(true_xy))
            continue;
        end
        if isempty(lr.est_pos_xy) || any(~isfinite(lr.est_pos_xy))
            continue;
        end

        n_probe = n_probe + 1;
        other_idx = find_matching_loc_result(other, lr.event_id);
        improved = false;
        if other_idx > 0
            lr_other = other(other_idx);
            err_primary = norm(lr.est_pos_xy - true_xy);
            err_other = norm(lr_other.est_pos_xy - true_xy);
            if strcmpi(mode_name, 'after')
                improved = err_primary < err_other;
            else
                improved = err_other < err_primary;
            end
            if improved
                n_improved = n_improved + 1;
            end
        end

        if strcmpi(mode_name, 'after')
            if improved
                est_color = [0.0 0.65 0.20];
                true_color = [0.0 0.55 0.15];
                line_color = [0.25 0.75 0.35];
            else
                est_color = [0.05 0.45 0.9];
                true_color = [0.1 0.1 0.1];
                line_color = [0.45 0.7 0.95];
            end
            est_marker = '+';
        else
            est_color = [0.75 0.2 0.2];
            true_color = [0.1 0.1 0.1];
            line_color = [0.9 0.55 0.55];
            est_marker = 'x';
        end

        plot(true_xy(1), true_xy(2), 'o', 'Color', true_color, ...
            'MarkerFaceColor', true_color, 'MarkerSize', 4);
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), est_marker, ...
            'Color', est_color, 'MarkerSize', 7, 'LineWidth', 1.3);
        plot([true_xy(1), lr.est_pos_xy(1)], [true_xy(2), lr.est_pos_xy(2)], ...
            '-', 'Color', line_color, 'LineWidth', 0.7);
    end
end


function idx = find_matching_loc_result(results, event_id)
    idx = 0;
    if isempty(results), return; end
    ids = [results.event_id];
    hit = find(ids == event_id, 1);
    if ~isempty(hit), idx = hit; end
end


function plot_target_panel(LocResults, events, GridValid, APs, ttl)
    plot(GridValid.xy(:,1), GridValid.xy(:,2), '.', 'Color', [0.94 0.94 0.94], 'MarkerSize', 2);
    hold on;
    plot(APs.pos_xy(:,1), APs.pos_xy(:,2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
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
        end
    end
    overlay_opportunistic_priors(events, true);
    axis equal tight;
    grid on;
    title(ttl);
    xlabel('X (m)');
    ylabel('Y (m)');
    add_cd_target_legend();
end


function plot_band_loc_results(LocResults, band_id)
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
            'Color', [0.9 0.15 0.15], 'MarkerFaceColor', [0.9 0.15 0.15], 'MarkerSize', 5);
        plot(lr.est_pos_xy(1), lr.est_pos_xy(2), 'x', ...
            'Color', [0.55 0 0], 'MarkerSize', 7, 'LineWidth', 1.3);
        plot([lr.true_pos_xy(1), lr.est_pos_xy(1)], ...
             [lr.true_pos_xy(2), lr.est_pos_xy(2)], ...
            '-', 'Color', [0.85 0.45 0.45], 'LineWidth', 0.8);
    end
end


function h = overlay_opportunistic_priors(events, with_label)
    h = gobjects(0);
    if isempty(events), return; end
    for i = 1:numel(events)
        lp = events(i).location_prior;
        b = events(i).band_id;
        switch lower(lp.type)
            case 'region'
                bbox = prior_bbox(lp);
                h(end+1) = rectangle('Position', [bbox(1), bbox(3), bbox(2)-bbox(1), bbox(4)-bbox(3)], ...
                    'EdgeColor', [0 0.45 0.2], 'LineWidth', 1.1); %#ok<AGROW>
                xy = [0.5*(bbox(1)+bbox(2)), 0.5*(bbox(3)+bbox(4))];
                plot(xy(1), xy(2), 's', 'Color', [0 0.5 0.2], 'MarkerSize', 5);
                tag = 'region';
            case 'gaussian'
                [mu, Sigma] = prior_gaussian(lp);
                h(end+1) = plot_gaussian_ellipse(mu, Sigma, [0.05 0.45 0.85]); %#ok<AGROW>
                plot(mu(1), mu(2), 'o', 'Color', [0.05 0.25 0.7], ...
                    'MarkerFaceColor', [0.05 0.25 0.7], 'MarkerSize', 4);
                xy = mu;
                tag = 'gauss';
            otherwise
                xy = representative_prior_xy(lp);
                tag = lp.type;
        end
        if with_label && all(isfinite(xy))
            text(xy(1) + 2, xy(2) + 2, sprintf('B%d CD %s', b, tag), ...
                'Color', [0 0.35 0.15], 'FontSize', 7, ...
                'BackgroundColor', [1 1 1 0.55]);
        end
    end
end


function h = plot_gaussian_ellipse(mu, Sigma, color)
    th = linspace(0, 2*pi, 120);
    [V, D] = eig(0.5 * (Sigma + Sigma'));
    A = V * sqrt(max(D, 0)) * sqrt(5.991);
    pts = mu(:) + A * [cos(th); sin(th)];
    h = plot(pts(1,:), pts(2,:), '--', 'Color', color, 'LineWidth', 1.0);
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


function xy = representative_prior_xy(lp)
    switch lower(lp.type)
        case 'region'
            bbox = prior_bbox(lp);
            xy = [0.5 * (bbox(1) + bbox(2)), 0.5 * (bbox(3) + bbox(4))];
        case 'gaussian'
            [mu, ~] = prior_gaussian(lp);
            xy = mu;
        case 'exact'
            if isfield(lp, 'xy'), xy = lp.xy(:)';
            elseif isfield(lp, 'value'), xy = lp.value(:)';
            else, xy = [NaN NaN];
            end
        otherwise
            xy = [NaN NaN];
    end
end


function bbox = prior_bbox(lp)
    if isfield(lp, 'bbox')
        bbox = lp.bbox(:)';
    elseif isfield(lp, 'value') && isnumeric(lp.value)
        bbox = lp.value(:)';
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'bbox')
        bbox = lp.value.bbox(:)';
    else
        bbox = [NaN NaN NaN NaN];
    end
    bbox = [min(bbox(1), bbox(2)), max(bbox(1), bbox(2)), ...
        min(bbox(3), bbox(4)), max(bbox(3), bbox(4))];
end


function [mu, Sigma] = prior_gaussian(lp)
    if isfield(lp, 'mu')
        mu = lp.mu(:)';
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'mu')
        mu = lp.value.mu(:)';
    else
        mu = [NaN NaN];
    end
    if isfield(lp, 'Sigma')
        Sigma = lp.Sigma;
    elseif isfield(lp, 'sigma')
        Sigma = eye(2) * lp.sigma^2;
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'Sigma')
        Sigma = lp.value.Sigma;
    elseif isfield(lp, 'value') && isstruct(lp.value) && isfield(lp.value, 'sigma')
        Sigma = eye(2) * lp.value.sigma^2;
    else
        Sigma = eye(2);
    end
    if isscalar(Sigma)
        Sigma = eye(2) * Sigma;
    elseif isvector(Sigma)
        Sigma = diag(Sigma(:));
    end
    Sigma = 0.5 * (Sigma + Sigma');
end


function F = get_band_dbm(band_fp)
    if isfield(band_fp, 'F_dBm')
        F = band_fp.F_dBm;
    elseif isfield(band_fp, 'RF_raw')
        F = band_fp.RF_raw;
    else
        error('[CD-demo] band fingerprint missing F_dBm/RF_raw');
    end
end


function add_cd_target_legend()
    h1 = plot(NaN, NaN, 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    h2 = plot(NaN, NaN, 'o', 'Color', [0.9 0.15 0.15], ...
        'MarkerFaceColor', [0.9 0.15 0.15], 'MarkerSize', 5);
    h3 = plot(NaN, NaN, 'x', 'Color', [0.55 0 0], 'MarkerSize', 7, 'LineWidth', 1.3);
    h4 = plot(NaN, NaN, 's', 'Color', [0 0.5 0.2], 'MarkerSize', 5);
    h5 = plot(NaN, NaN, 'o', 'Color', [0.05 0.25 0.7], ...
        'MarkerFaceColor', [0.05 0.25 0.7], 'MarkerSize', 4);
    legend([h1 h2 h3 h4 h5], {'AP','target true','target est','CD region opp','CD gaussian opp'}, ...
        'Location', 'best');
end


function add_cd_probe_legend(mode_name)
    h1 = plot(NaN, NaN, 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    h2 = plot(NaN, NaN, 'o', 'Color', [0.1 0.1 0.1], ...
        'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 4);
    h3 = plot(NaN, NaN, 's', 'Color', [0 0.5 0.2], 'MarkerSize', 5);
    h4 = plot(NaN, NaN, 'o', 'Color', [0.05 0.25 0.7], ...
        'MarkerFaceColor', [0.05 0.25 0.7], 'MarkerSize', 4);
    if strcmpi(mode_name, 'after')
        h5 = plot(NaN, NaN, '+', 'Color', [0.0 0.65 0.20], 'MarkerSize', 7, 'LineWidth', 1.3);
        h6 = plot(NaN, NaN, '+', 'Color', [0.05 0.45 0.9], 'MarkerSize', 7, 'LineWidth', 1.3);
        legend([h1 h2 h5 h6 h3 h4], ...
            {'AP','probe true','after improved','after not improved','CD region opp','CD gaussian opp'}, ...
            'Location', 'best');
    else
        h5 = plot(NaN, NaN, 'x', 'Color', [0.75 0.2 0.2], 'MarkerSize', 7, 'LineWidth', 1.2);
        legend([h1 h2 h5 h3 h4], ...
            {'AP','probe true','before est','CD region opp','CD gaussian opp'}, ...
            'Location', 'best');
    end
end
