function [M0State, BandWinners] = resolve_band_occupancy_m0( ...
    BroadbandCandidates, OppCandidatesPerBand, TargetCandidatesPerBand, ...
    PersistentCandidates, M0State, Config)
% RESOLVE_BAND_OCCUPANCY_M0  统一优先级裁决（信源发生器内部完成）
%
%   优先级: broadband_cal(4) > opportunistic(3) > target(2) > persistent_cal(1)
%   后续模块不得再做优先级判定。

    B = Config.m0.num_bands;

    priority_of = struct( ...
        'broadband_cal',   4, ...
        'opportunistic',   3, ...
        'target',          2, ...
        'persistent_cal',  1);

    BandWinners = struct('has_source', cell(1, B), 'winner', cell(1, B));

    for b = 1:B
        all_cands = [];

        % 1) broadband_cal（覆盖全频带）
        for j = 1:numel(BroadbandCandidates)
            bc = BroadbandCandidates(j);
            if any(bc.band_list == b)
                bc.priority = priority_of.broadband_cal;
                all_cands = append_cand(all_cands, normalize_cand(bc));
            end
        end

        % 2) opportunistic
        if ~isempty(OppCandidatesPerBand(b).candidates)
            for ci = 1:numel(OppCandidatesPerBand(b).candidates)
                oc = OppCandidatesPerBand(b).candidates(ci);
                oc.priority = priority_of.opportunistic;
                all_cands = append_cand(all_cands, normalize_cand(oc));
            end
        end

        % 3) target
        if ~isempty(TargetCandidatesPerBand(b).candidates)
            for ci = 1:numel(TargetCandidatesPerBand(b).candidates)
                tc = TargetCandidatesPerBand(b).candidates(ci);
                tc.priority = priority_of.target;
                all_cands = append_cand(all_cands, normalize_cand(tc));
            end
        end

        % 4) persistent_cal（保底）
        for j = 1:numel(PersistentCandidates)
            pc = PersistentCandidates(j);
            if any(pc.band_list == b) || pc.band_id == b
                pc.priority = priority_of.persistent_cal;
                all_cands = append_cand(all_cands, normalize_cand(pc));
            end
        end

        if isempty(all_cands)
            BandWinners(b).has_source = false;
            BandWinners(b).winner     = [];
            M0State.active_per_band(b).has_source   = false;
            M0State.active_per_band(b).instance_id  = 0;
            M0State.active_per_band(b).template_key = '';
        else
            N = numel(all_cands);
            keys = zeros(N, 4);
            for ci = 1:N
                keys(ci, 1) = all_cands(ci).priority;
                keys(ci, 2) = double(all_cands(ci).is_continuing);
                keys(ci, 3) = all_cands(ci).true_tx_power;
                keys(ci, 4) = -all_cands(ci).template_idx;
            end
            [~, order] = sortrows(keys, [-1, -2, -3, -4]);
            winner = all_cands(order(1));
            winner.is_selected_final = true;

            BandWinners(b).has_source = true;
            BandWinners(b).winner     = winner;
            M0State.active_per_band(b).has_source   = true;
            M0State.active_per_band(b).instance_id  = winner.instance_id;
            M0State.active_per_band(b).template_key = winner.template_key;
        end
    end
end


function c = normalize_cand(c)
% 确保所有候选具有统一字段集合
    defaults = {
        'instance_id',      0; ...
        'template_key',     ''; ...
        'source_type',      ''; ...
        'source_subtype',   ''; ...
        'band_list',        []; ...
        'band_id',          0; ...
        'true_position',    [0 0]; ...
        'true_tx_power',    0; ...
        'tx_power_by_band', []; ...
        'is_continuing',    false; ...
        'template_idx',     0; ...
        'is_persistent',    false; ...
        'priority',         0; ...
        'location_prior',   struct('type','none','value',[]); ...
        'power_prior',      struct('type','none','value',[]) ...
    };
    for i = 1:size(defaults, 1)
        fn = defaults{i, 1};
        if ~isfield(c, fn)
            c.(fn) = defaults{i, 2};
        end
    end
end


function arr = append_cand(arr, c)
    if isempty(arr)
        arr = c;
    else
        arr(end+1) = c;
    end
end
