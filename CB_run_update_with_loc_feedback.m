function [APProfile_best, SpatialFP_best, diagnostics_best] = CB_run_update_with_loc_feedback( ...
    APProfile_current, SpatialFP_current, calibration_events, validation_events, FrameStates, Config)
% CB_RUN_UPDATE_WITH_LOC_FEEDBACK  CB update with optional Loc-RMSE GA feedback.

    if nargin < 6
        Config = struct();
    end
    cfg = CB_calib_ap_profile_defaults(Config);

    [APProfile_base, SpatialFP_base, diag_base] = CB_update_ap_profiles_with_calibration_sources( ...
        APProfile_current, SpatialFP_current, calibration_events, Config);

    diagnostics_best = diag_base;
    APProfile_best = APProfile_base;
    SpatialFP_best = SpatialFP_base;
    diagnostics_best.feedback = struct('enabled', false, 'reason', 'not_requested');

    if ~cfg.enable_ga_feedback
        diagnostics_best.feedback.reason = 'disabled';
        return;
    end
    if nargin < 4 || isempty(validation_events) || nargin < 5 || isempty(FrameStates)
        diagnostics_best.feedback.reason = 'no_validation_events_or_framestates';
        return;
    end

    rmse0 = compute_loc_rmse(validation_events, SpatialFP_base, FrameStates, Config);
    diagnostics_best.feedback.enabled = true;
    diagnostics_best.feedback.initial_loc_rmse_m = rmse0;
    diagnostics_best.feedback.threshold_m = cfg.loc_rmse_threshold_m;

    if ~isfinite(rmse0) || rmse0 <= cfg.loc_rmse_threshold_m
        diagnostics_best.feedback.reason = 'initial_rmse_accepted';
        return;
    end

    [cfg_best, rmse_best, hist] = run_ga_search(cfg, APProfile_current, SpatialFP_current, ...
        calibration_events, validation_events, FrameStates, Config);

    Config_best = Config;
    Config_best.calib_ap_profile = cfg_best;
    [APProfile_best, SpatialFP_best, diagnostics_best] = CB_update_ap_profiles_with_calibration_sources( ...
        APProfile_current, SpatialFP_current, calibration_events, Config_best);
    diagnostics_best.feedback.enabled = true;
    diagnostics_best.feedback.reason = 'ga_feedback_applied';
    diagnostics_best.feedback.initial_loc_rmse_m = rmse0;
    diagnostics_best.feedback.best_loc_rmse_m = rmse_best;
    diagnostics_best.feedback.ga_history = hist;
    diagnostics_best.feedback.best_config = cfg_best;
end


function [cfg_best, rmse_best, hist] = run_ga_search(cfg, APProfile_current, SpatialFP_current, ...
    calibration_events, validation_events, FrameStates, Config)

    pop_n = max(4, cfg.ga_population);
    ngen = max(1, cfg.ga_generations);
    elite_n = min(max(1, cfg.ga_elite_count), pop_n);
    genes0 = cfg_to_gene(cfg);
    n_gene = numel(genes0);

    pop = zeros(pop_n, n_gene);
    pop(1, :) = genes0;
    for i = 2:pop_n
        pop(i, :) = mutate_gene(genes0, cfg.ga_mutation_rate, 0.35);
    end

    hist = struct('best_rmse', zeros(1, ngen), 'best_gene', zeros(ngen, n_gene));
    rmse = inf(pop_n, 1);
    for gen = 1:ngen
        for i = 1:pop_n
            cfg_i = gene_to_cfg(cfg, pop(i, :));
            Config_i = Config;
            Config_i.calib_ap_profile = cfg_i;
            try
                [~, SpatialFP_i] = CB_update_ap_profiles_with_calibration_sources( ...
                    APProfile_current, SpatialFP_current, calibration_events, Config_i);
                rmse(i) = compute_loc_rmse(validation_events, SpatialFP_i, FrameStates, Config_i);
            catch ME
                fprintf('[CB-GA] candidate failed: %s\n', ME.message);
                rmse(i) = inf;
            end
        end

        [rmse_sorted, ord] = sort(rmse, 'ascend');
        pop = pop(ord, :);
        rmse = rmse_sorted;
        hist.best_rmse(gen) = rmse(1);
        hist.best_gene(gen, :) = pop(1, :);

        new_pop = pop;
        for i = (elite_n + 1):pop_n
            p1 = tournament(pop, rmse);
            p2 = tournament(pop, rmse);
            child = crossover_gene(p1, p2, cfg.ga_crossover_rate);
            child = mutate_gene(child, cfg.ga_mutation_rate, 0.25);
            new_pop(i, :) = child;
        end
        pop = new_pop;
    end

    cfg_best = gene_to_cfg(cfg, hist.best_gene(end, :));
    rmse_best = hist.best_rmse(end);
end


function rmse = compute_loc_rmse(events, SpatialFP, FrameStates, Config)
    ev = events;
    for i = 1:numel(ev)
        ev(i).route_action = 'localize_only';
        if ~isfield(ev(i), 'type_hat') || isempty(ev(i).type_hat)
            ev(i).type_hat = 'cb_validation';
        end
    end
    LocResults = run_m4_wknn_localization(ev, SpatialFP, FrameStates, Config);
    if isempty(LocResults) || ~isfield(LocResults, 'loc_error')
        rmse = inf;
        return;
    end
    err = [LocResults.loc_error];
    err = err(isfinite(err));
    if isempty(err)
        rmse = inf;
    else
        rmse = sqrt(mean(err.^2));
    end
end


function g = cfg_to_gene(cfg)
    g = log([cfg.sigma_f_dB, cfg.ell_x_m, cfg.ell_y_m, cfg.sigma0_dB, cfg.sigma_meas_dB]);
end


function cfg2 = gene_to_cfg(cfg, g)
    vals = exp(g);
    cfg2 = cfg;
    cfg2.sigma_f_dB = clamp(vals(1), 0.5, 30);
    cfg2.ell_x_m = clamp(vals(2), 3, 250);
    cfg2.ell_y_m = clamp(vals(3), 3, 250);
    cfg2.sigma0_dB = clamp(vals(4), 0.5, 30);
    cfg2.sigma_meas_dB = clamp(vals(5), 0.1, 15);
end


function child = crossover_gene(p1, p2, rate)
    if rand() >= rate
        child = p1;
        return;
    end
    mask = rand(size(p1)) < 0.5;
    child = p1;
    child(mask) = p2(mask);
end


function g2 = mutate_gene(g, rate, scale)
    g2 = g;
    mask = rand(size(g)) < rate;
    g2(mask) = g2(mask) + scale * randn(size(g2(mask)));
end


function p = tournament(pop, rmse)
    n = size(pop, 1);
    a = randi(n);
    b = randi(n);
    if rmse(a) <= rmse(b)
        p = pop(a, :);
    else
        p = pop(b, :);
    end
end


function y = clamp(x, lo, hi)
    y = min(max(x, lo), hi);
end
