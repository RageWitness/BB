function [template_scores, best_template_key, best_score, second_best_score, diag_info] = ...
    compute_group_position_match_m3(GroupItem, EventList, SpatialFP, TemplateList, Config)
% COMPUTE_GROUP_POSITION_MATCH_M3  group 级逐模板位置匹配
%
%   S_pos(j) = max_{g in P_j} sum_{b in B_obs} [omega_b * exp(-(D_b(g)^2)/(tau_b^2))]

    B = numel(SpatialFP.band);
    G = size(SpatialFP.grid_xy, 1);
    cfg = fill_pos_match_defaults(Config, B);

    template_scores = zeros(1, numel(TemplateList));
    best_template_key = '';
    best_score = 0;
    second_best_score = 0;

    diag_info = struct();
    diag_info.band_distance_to_grid = nan(B, G);
    diag_info.band_best_grid_idx = nan(1, B);
    diag_info.band_best_pos_xy = nan(B, 2);
    diag_info.multiband_position_consistency = NaN;
    diag_info.omega_per_band = cfg.omega_per_band;
    diag_info.tau_per_band = cfg.tau_per_band;

    if isempty(EventList) || isempty(get_field_default(GroupItem, 'active_bands', []))
        return;
    end

    active_bands = GroupItem.active_bands(:).';

    % 1) 预计算 D_b(g)
    for b = active_bands
        eidx = GroupItem.member_event_idx_per_band(b);
        if eidx <= 0 || eidx > numel(EventList)
            continue;
        end
        ev = EventList(eidx);
        [~, s_obs] = aggregate_event_fingerprint_m4(ev);
        F_shape = SpatialFP.band(b).F_shape_l1;
        dvec = sqrt(sum((s_obs - F_shape).^2, 1));
        diag_info.band_distance_to_grid(b, :) = dvec;
        [~, g_best] = min(dvec);
        diag_info.band_best_grid_idx(b) = g_best;
        diag_info.band_best_pos_xy(b, :) = SpatialFP.grid_xy(g_best, :);
    end

    % 2) 多 band 位置离散度
    valid = ~isnan(diag_info.band_best_grid_idx(active_bands));
    if sum(valid) >= 2
        pos_use = diag_info.band_best_pos_xy(active_bands(valid), :);
        c = mean(pos_use, 1);
        d = sqrt(sum((pos_use - c).^2, 2));
        diag_info.multiband_position_consistency = max(d);
    end

    if isempty(TemplateList)
        return;
    end

    % 3) 逐模板评分
    for j = 1:numel(TemplateList)
        tpl = TemplateList(j);
        cand_pos = get_candidate_positions(tpl);
        if isempty(cand_pos) || size(cand_pos, 2) ~= 2
            template_scores(j) = 0;
            continue;
        end

        [cand_grid_idx, ~] = map_positions_to_grid(cand_pos, SpatialFP.grid_xy);
        cand_grid_idx = unique(cand_grid_idx(:).');
        cand_grid_idx = cand_grid_idx(cand_grid_idx >= 1 & cand_grid_idx <= G);
        if isempty(cand_grid_idx)
            template_scores(j) = 0;
            continue;
        end

        band_use = select_bands_for_template(tpl, active_bands, B);
        if isempty(band_use)
            template_scores(j) = 0;
            continue;
        end

        s_best = 0;
        wsum = sum(cfg.omega_per_band(band_use)) + cfg.eps_val;
        for gg = cand_grid_idx
            s = 0;
            for bb = band_use
                d_bg = diag_info.band_distance_to_grid(bb, gg);
                if isnan(d_bg) || isinf(d_bg)
                    continue;
                end
                s = s + cfg.omega_per_band(bb) * exp(-(d_bg.^2) / (cfg.tau_per_band(bb).^2 + cfg.eps_val));
            end
            s_best = max(s_best, s / wsum);
        end
        template_scores(j) = min(max(s_best, 0), 1);
    end

    [best_score, idx_best] = max(template_scores);
    if ~isempty(idx_best) && idx_best >= 1 && idx_best <= numel(TemplateList)
        best_template_key = get_field_default(TemplateList(idx_best), 'template_key', '');
    end

    if numel(template_scores) >= 2
        ss = sort(template_scores, 'descend');
        second_best_score = ss(2);
    else
        second_best_score = 0;
    end
end


%% ==================== 局部函数 ====================

function cfg = fill_pos_match_defaults(Config, B)
    cfg.tau_per_band = 0.35 * ones(1, B);
    cfg.omega_per_band = ones(1, B);
    cfg.eps_val = 1e-12;

    if isfield(Config, 'm3') && isfield(Config.m3, 'position_match')
        pm = Config.m3.position_match;
        if isfield(pm, 'tau_per_band')
            tau = pm.tau_per_band;
            if isscalar(tau), tau = repmat(tau, 1, B); end
            if numel(tau) >= B, cfg.tau_per_band = tau(1:B); end
        end
        if isfield(pm, 'omega_per_band')
            om = pm.omega_per_band;
            if isscalar(om), om = repmat(om, 1, B); end
            if numel(om) >= B, cfg.omega_per_band = om(1:B); end
        end
    end
end


function [idx, dist] = map_positions_to_grid(pos_xy, grid_xy)
    K = size(pos_xy, 1);
    idx = zeros(K, 1);
    dist = zeros(K, 1);
    for k = 1:K
        d = sqrt(sum((grid_xy - pos_xy(k, :)).^2, 2));
        [dist(k), idx(k)] = min(d);
    end
end


function cand_pos = get_candidate_positions(tpl)
    cand_pos = [];
    if isfield(tpl, 'position_prior') && isstruct(tpl.position_prior)
        cand_pos = get_field_default(tpl.position_prior, 'candidate_positions', []);
    end
    if isempty(cand_pos)
        cand_pos = get_field_default(tpl, 'candidate_positions', []);
    end
    if isempty(cand_pos) && isfield(tpl, 'fixed_pos_xy')
        cand_pos = tpl.fixed_pos_xy;
    end
end


function bands = select_bands_for_template(tpl, active_bands, B)
    bands = active_bands;
    bm = get_field_default(tpl, 'band_mask', []);
    if isempty(bm)
        return;
    end
    bm = bm(:).';
    if numel(bm) < B
        bm = [bm, zeros(1, B - numel(bm))];
    end
    bands = active_bands(bm(active_bands) > 0);
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
