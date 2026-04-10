function [M0State, BandWinners] = resolve_band_occupancy_m0( ...
    TrustedCandidates, PriorCandidatesPerBand, TargetCandidatesPerBand, ...
    M0State, Config)
% RESOLVE_BAND_OCCUPANCY_M0  对每个频带裁决最终唯一主导源
%
%   候选和 winner 均携带完整先验对象（position_prior, time_prior 等）。
%   最终保证：N_active(b,t) <= 1

    B = Config.m0.num_bands;

    priority_map = containers.Map( ...
        {'trusted_fixed', 'prior_pos_known', 'prior_time_known', 'ordinary_target'}, ...
        {4, 3, 2, 1});

    for b = 1:B
        all_candidates = [];

        % --- 1. 收集 trusted 候选 ---
        for j = 1:numel(TrustedCandidates)
            tc = TrustedCandidates(j);
            if tc.is_active && ismember(b, tc.bands_covered)
                c = pack_candidate(tc, 'trusted_fixed', 'trusted_fixed', ...
                    priority_map('trusted_fixed'), true);
                all_candidates = append_cand(all_candidates, c);
            end
        end

        % --- 2. 收集 prior 候选 ---
        if ~isempty(PriorCandidatesPerBand(b).candidates)
            for ci = 1:numel(PriorCandidatesPerBand(b).candidates)
                pc = PriorCandidatesPerBand(b).candidates(ci);
                c = pack_candidate(pc, pc.source_type, pc.source_subtype, ...
                    priority_map(pc.source_subtype), pc.is_continuing);
                all_candidates = append_cand(all_candidates, c);
            end
        end

        % --- 3. 收集 target 候选 ---
        if ~isempty(TargetCandidatesPerBand(b).candidates)
            for ci = 1:numel(TargetCandidatesPerBand(b).candidates)
                tc = TargetCandidatesPerBand(b).candidates(ci);
                c = pack_candidate(tc, 'ordinary_target', 'ordinary_target', ...
                    priority_map('ordinary_target'), tc.is_continuing);
                all_candidates = append_cand(all_candidates, c);
            end
        end

        % --- 4. 裁决 ---
        if isempty(all_candidates)
            BandWinners(b).has_source = false;
            BandWinners(b).winner     = [];
            M0State.active_per_band(b).has_source   = false;
            M0State.active_per_band(b).instance_id  = 0;
            M0State.active_per_band(b).template_key = '';
        else
            N_cand = numel(all_candidates);
            sort_keys = zeros(N_cand, 4);
            for ci = 1:N_cand
                sort_keys(ci, 1) = all_candidates(ci).priority;
                sort_keys(ci, 2) = double(all_candidates(ci).is_continuing);
                sort_keys(ci, 3) = all_candidates(ci).tx_power_dBm;
                sort_keys(ci, 4) = -all_candidates(ci).template_idx;
            end
            [~, order] = sortrows(sort_keys, [-1, -2, -3, -4]);
            winner = all_candidates(order(1));

            BandWinners(b).has_source = true;
            BandWinners(b).winner     = winner;

            M0State.active_per_band(b).has_source   = true;
            M0State.active_per_band(b).instance_id  = winner.instance_id;
            M0State.active_per_band(b).template_key = winner.template_key;
        end
    end
end


function c = pack_candidate(src, source_type, source_subtype, priority, is_continuing)
% PACK_CANDIDATE  将候选源打包为统一格式（含先验对象）
    c.instance_id           = src.instance_id;
    c.template_key          = src.template_key;
    c.source_type           = source_type;
    c.source_subtype        = source_subtype;
    c.true_pos_xy           = src.true_pos_xy;
    c.tx_power_dBm          = src.tx_power_dBm;
    c.priority              = priority;
    c.is_continuing         = is_continuing;
    c.template_idx          = src.template_idx;
    c.position_prior        = src.position_prior;
    c.time_prior            = src.time_prior;
    c.power_nominal_dBm     = src.power_nominal_dBm;
    c.power_range_dBm       = src.power_range_dBm;
    c.power_stability_level = src.power_stability_level;
    c.credibility_prior     = src.credibility_prior;
    c.upgrade_potential     = src.upgrade_potential;
end


function arr = append_cand(arr, c)
% APPEND_CAND  向候选数组追加一个元素
    if isempty(arr)
        arr = c;
    else
        arr(end+1) = c; %#ok<AGROW>
    end
end
