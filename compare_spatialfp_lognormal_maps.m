%% COMPARE_SPATIALFP_LOGNORMAL_MAPS
% Compare two lognormal SpatialFP libraries as GPR residual targets.
%
% Main view:
%   Delta_b,m(g) = SpatialFP_new.band(b).RF_raw(m,g)
%                - SpatialFP_old.band(b).RF_raw(m,g)
%
% This is exactly the map residual that CA/GPR would try to learn when the
% online map is the old library and the new library is treated as reference.

clear; clc; close all;

old_file = fullfile('cache', 'SpatialFP_lognormal.mat');
new_file = fullfile('cache', 'SpatialFP_lognormal_new.mat');

gpr_field = 'RF_raw';          % CA now calibrates RF_raw, i.e. raw RSS dBm
show_ap_detail = true;         % per-band AP residual heatmaps
save_figures = false;
out_dir = fullfile('figures', 'fp_compare_lognormal');

if save_figures && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('============================================\n');
fprintf('  SpatialFP lognormal library comparison\n');
fprintf('============================================\n\n');
fprintf('Old: %s\n', old_file);
fprintf('New: %s\n\n', new_file);

SpatialFP_old = load_spatialfp(old_file);
SpatialFP_new = load_spatialfp(new_file);

check_compatible(SpatialFP_old, SpatialFP_new);

B = SpatialFP_old.B;
M = SpatialFP_old.M;
G = SpatialFP_old.G;
grid_xy = SpatialFP_old.grid_xy;

fprintf('B=%d bands, M=%d AP, G=%d grid points\n\n', B, M, G);

%% 1. Field-level summary
fields = {'RF_raw', 'F_dBm', 'centered_dBm', 'RF_minmax', 'F_shape_l1'};
fprintf('--- Field delta summary: new - old ---\n');
fprintf('%-14s %-6s %-12s %-12s %-12s %-12s %-12s %-10s\n', ...
    'Field', 'Band', 'meanAbs', 'rms', 'p95Abs', 'maxAbs', 'meanDelta', 'corr');

Summary = struct([]);
si = 0;
for fi = 1:numel(fields)
    fld = fields{fi};
    for b = 1:B
        if ~isfield(SpatialFP_old.band(b), fld) || ~isfield(SpatialFP_new.band(b), fld)
            continue;
        end
        A = SpatialFP_old.band(b).(fld);
        C = SpatialFP_new.band(b).(fld);
        D = C - A;
        vals = D(:);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end

        corr_val = corr_safe(A(:), C(:));
        si = si + 1;
        Summary(si).field = fld; %#ok<SAGROW>
        Summary(si).band = b;
        Summary(si).mean_abs = mean(abs(vals));
        Summary(si).rms = sqrt(mean(vals.^2));
        Summary(si).p95_abs = prctile(abs(vals), 95);
        Summary(si).max_abs = max(abs(vals));
        Summary(si).mean_delta = mean(vals);
        Summary(si).corr = corr_val;

        fprintf('%-14s %-6d %-12.4g %-12.4g %-12.4g %-12.4g %-12.4g %-10.4f\n', ...
            fld, b, Summary(si).mean_abs, Summary(si).rms, Summary(si).p95_abs, ...
            Summary(si).max_abs, Summary(si).mean_delta, Summary(si).corr);
    end
end
fprintf('\n');

%% 2. GPR target view: Delta RF_raw, averaged over AP
DeltaBand = cell(1, B);
MeanDelta = cell(1, B);
RmsDelta = cell(1, B);
all_mean = [];
all_rms = [];

for b = 1:B
    D = get_field(SpatialFP_new.band(b), gpr_field) - ...
        get_field(SpatialFP_old.band(b), gpr_field);
    DeltaBand{b} = D;
    MeanDelta{b} = mean(D, 1);
    RmsDelta{b} = sqrt(mean(D.^2, 1));
    all_mean = [all_mean; MeanDelta{b}(:)]; %#ok<AGROW>
    all_rms = [all_rms; RmsDelta{b}(:)]; %#ok<AGROW>
end

mean_clim = symmetric_clim(all_mean);
rms_clim = [0, max(all_rms(isfinite(all_rms)))];
if ~isfinite(rms_clim(2)) || rms_clim(2) <= 0, rms_clim = [0, 1]; end

figure('Name', 'GPR target mean delta RF_raw', 'Position', [50 60 1350 380]);
for b = 1:B
    subplot(1, B, b);
    scatter(grid_xy(:,1), grid_xy(:,2), 12, MeanDelta{b}, 'filled');
    caxis(mean_clim); colorbar;
    axis equal tight; grid on;
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('Band %d mean_m Delta RF\\_raw', b));
end
sgtitle('GPR target: mean over AP of Delta RF\_raw = new - old (dB)');
save_current_fig(save_figures, out_dir, '01_mean_delta_rf_raw');

figure('Name', 'GPR target RMS delta RF_raw', 'Position', [60 90 1350 380]);
for b = 1:B
    subplot(1, B, b);
    scatter(grid_xy(:,1), grid_xy(:,2), 12, RmsDelta{b}, 'filled');
    caxis(rms_clim); colorbar;
    axis equal tight; grid on;
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('Band %d RMS_m Delta RF\\_raw', b));
end
sgtitle('GPR target: AP RMS of Delta RF\_raw (dB)');
save_current_fig(save_figures, out_dir, '02_rms_delta_rf_raw');

%% 3. Per-AP residual strength, which equals GPR output columns Y(:,m)
mean_abs_ap = zeros(B, M);
rms_ap = zeros(B, M);
for b = 1:B
    D = DeltaBand{b};
    mean_abs_ap(b, :) = mean(abs(D), 2)';
    rms_ap(b, :) = sqrt(mean(D.^2, 2))';
end

figure('Name', 'GPR Y columns per AP strength', 'Position', [100 120 1100 420]);
subplot(1,2,1);
bar(mean_abs_ap');
grid on; xlabel('AP index'); ylabel('mean |Delta RF\_raw| (dB)');
title('Mean absolute GPR target per AP');
legend(make_band_labels(B), 'Location', 'best');

subplot(1,2,2);
bar(rms_ap');
grid on; xlabel('AP index'); ylabel('RMS Delta RF\_raw (dB)');
title('RMS GPR target per AP');
legend(make_band_labels(B), 'Location', 'best');
save_current_fig(save_figures, out_dir, '03_per_ap_delta_strength');

%% 4. Histograms and old-vs-new scatter
figure('Name', 'Delta RF_raw histogram and old-vs-new', 'Position', [120 140 1300 650]);
for b = 1:B
    D = DeltaBand{b};
    oldF = get_field(SpatialFP_old.band(b), gpr_field);
    newF = get_field(SpatialFP_new.band(b), gpr_field);

    subplot(2, B, b);
    histogram(D(:), 50);
    grid on;
    xlabel('Delta RF\_raw (dB)'); ylabel('count');
    title(sprintf('Band %d delta histogram', b));

    subplot(2, B, B + b);
    sample_idx = sample_indices(numel(oldF), min(12000, numel(oldF)));
    scatter(oldF(sample_idx), newF(sample_idx), 5, 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    hold on;
    lims = axis_equal_limits(oldF(sample_idx), newF(sample_idx));
    plot(lims, lims, 'r--', 'LineWidth', 1);
    hold off;
    xlim(lims); ylim(lims); grid on; axis square;
    xlabel('old RF\_raw (dBm)'); ylabel('new RF\_raw (dBm)');
    title(sprintf('Band %d old vs new', b));
end
sgtitle('Distribution view of RF\_raw library change');
save_current_fig(save_figures, out_dir, '04_hist_old_vs_new');

%% 5. Smoothness view for GPR: spatial distance vs residual difference
figure('Name', 'Residual field spatial smoothness', 'Position', [140 160 1250 360]);
for b = 1:B
    subplot(1, B, b);
    Dmean = MeanDelta{b}(:);
    [d_space, d_delta] = sampled_pair_differences(grid_xy, Dmean, 20000);
    scatter(d_space, d_delta, 4, 'filled', ...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
    grid on;
    xlabel('spatial distance (m)');
    ylabel('|Delta(g_i)-Delta(g_j)| (dB)');
    title(sprintf('Band %d residual smoothness', b));
end
sgtitle('GPR difficulty view: smooth residual fields are easier to learn');
save_current_fig(save_figures, out_dir, '05_residual_smoothness');

%% 6. Optional per-band AP heatmaps
if show_ap_detail
    n_cols = ceil(sqrt(M));
    n_rows = ceil(M / n_cols);
    for b = 1:B
        D = DeltaBand{b};
        ap_clim = symmetric_clim(D(:));
        figure('Name', sprintf('Band %d per-AP Delta RF_raw', b), ...
            'Position', [80 + 20*b, 80 + 20*b, 260*n_cols, 230*n_rows]);
        for m = 1:M
            subplot(n_rows, n_cols, m);
            scatter(grid_xy(:,1), grid_xy(:,2), 10, D(m,:), 'filled');
            caxis(ap_clim); colorbar;
            axis equal tight; grid on;
            xlabel('X'); ylabel('Y');
            title(sprintf('B%d AP%d', b, m));
        end
        sgtitle(sprintf('Band %d: per-AP GPR target Delta RF\\_raw (dB)', b));
        save_current_fig(save_figures, out_dir, sprintf('06_band%d_per_ap_delta', b));
    end
end

fprintf('--- GPR training-form interpretation ---\n');
fprintf('For each band b and AP m, GPR would use:\n');
fprintf('  X = grid_xy                         size %d x 2\n', G);
fprintf('  Y(:,m) = Delta RF_raw(m,:)''         size %d x %d per band\n', G, M);
fprintf('A large/smooth Delta map is learnable by SE/Matérn kernels.\n');
fprintf('Sharp isolated spikes or sign-changing patches need local gating or shorter length scale.\n\n');

fprintf('Done. Figures generated in MATLAB. save_figures=%d\n', save_figures);


%% ==================== helpers ====================

function SpatialFP = load_spatialfp(fname)
    if ~exist(fname, 'file')
        error('File not found: %s', fname);
    end
    S = load(fname);
    if isfield(S, 'SpatialFP')
        SpatialFP = S.SpatialFP;
        return;
    end

    names = fieldnames(S);
    for i = 1:numel(names)
        v = S.(names{i});
        if isstruct(v) && isfield(v, 'band') && isfield(v, 'grid_xy')
            SpatialFP = v;
            return;
        end
    end
    error('No SpatialFP-like variable found in %s', fname);
end


function check_compatible(A, B)
    must_equal(A.B, B.B, 'B');
    must_equal(A.M, B.M, 'M');
    must_equal(A.G, B.G, 'G');
    if any(size(A.grid_xy) ~= size(B.grid_xy)) || max(abs(A.grid_xy(:) - B.grid_xy(:))) > 1e-9
        error('grid_xy differs between two libraries; compare requires same grid order.');
    end
end


function must_equal(a, b, name)
    if a ~= b
        error('SpatialFP.%s mismatch: old=%g, new=%g', name, a, b);
    end
end


function F = get_field(band_fp, fld)
    if isfield(band_fp, fld)
        F = band_fp.(fld);
    elseif strcmp(fld, 'RF_raw') && isfield(band_fp, 'F_dBm')
        F = band_fp.F_dBm;
    else
        error('Missing field %s in SpatialFP.band', fld);
    end
end


function c = corr_safe(a, b)
    a = a(:); b = b(:);
    ok = isfinite(a) & isfinite(b);
    if sum(ok) < 2
        c = NaN;
    else
        C = corrcoef(a(ok), b(ok));
        c = C(1,2);
    end
end


function clim = symmetric_clim(v)
    v = v(:);
    v = v(isfinite(v));
    if isempty(v)
        clim = [-1, 1];
        return;
    end
    a = max(abs(v));
    if a <= 0 || ~isfinite(a), a = 1; end
    clim = [-a, a];
end


function labels = make_band_labels(B)
    labels = cell(1, B);
    for b = 1:B
        labels{b} = sprintf('Band %d', b);
    end
end


function idx = sample_indices(N, K)
    if N <= K
        idx = 1:N;
    else
        idx = randperm(N, K);
    end
end


function lims = axis_equal_limits(x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        lims = [0, 1];
        return;
    end
    lo = min(vals); hi = max(vals);
    if lo == hi
        lo = lo - 0.5; hi = hi + 0.5;
    end
    pad = 0.05 * (hi - lo);
    lims = [lo - pad, hi + pad];
end


function [d_space, d_delta] = sampled_pair_differences(grid_xy, delta_vec, n_pair)
    G = size(grid_xy, 1);
    if G < 2
        d_space = [];
        d_delta = [];
        return;
    end
    n_pair = min(n_pair, G * max(G - 1, 1));
    i1 = randi(G, n_pair, 1);
    i2 = randi(G, n_pair, 1);
    same = i1 == i2;
    i2(same) = mod(i2(same), G) + 1;

    dxy = grid_xy(i1, :) - grid_xy(i2, :);
    d_space = sqrt(sum(dxy.^2, 2));
    d_delta = abs(delta_vec(i1) - delta_vec(i2));

    ok = isfinite(d_space) & isfinite(d_delta);
    d_space = d_space(ok);
    d_delta = d_delta(ok);
end


function save_current_fig(enable, out_dir, name)
    if ~enable, return; end
    saveas(gcf, fullfile(out_dir, [name, '.png']));
    savefig(gcf, fullfile(out_dir, [name, '.fig']));
end
