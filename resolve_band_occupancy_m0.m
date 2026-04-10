function [M0State, BandWinners] = resolve_band_occupancy_m0( ...
    TrustedCandidates, PriorCandidatesPerBand, TargetCandidatesPerBand, ...
    M0State, Config)
% RESOLVE_BAND_OCCUPANCY_M0  对每个频带裁决最终唯一主导源
%
%   [M0State, BandWinners] = resolve_band_occupancy_m0(...)
%
%   裁决优先级（从高到低）：
%     1. trusted_fixed
%     2. prior_pos_known
%     3. prior_time_known
%     4. ordinary_target
%
%   同类竞争规则：
%     1. 已在持续中的实例优先
%     2. 发射功率更高者优先
%     3. 模板编号更小者优先
%
%   最终保证：N_active(b,t) <= 1

    B = Config.m0.num_bands;

    % 优先级映射
    priority_map = containers.Map( ...
        {'trusted_fixed', 'prior_pos_known', 'prior_time_known', 'ordinary_target'}, ...
        {4, 3, 2, 1});

    for b = 1:B
        all_candidates = [];

        % --- 1. 收集 trusted 候选 ---
        for j = 1:numel(TrustedCandidates)
            tc = TrustedCandidates(j);
            if tc.is_active && ismember(b, tc.bands_covered)
                c.instance_id   = tc.instance_id;
                c.template_key  = tc.template_key;
                c.source_type   = 'trusted_fixed';
                c.source_subtype = 'trusted_fixed';
                c.true_pos_xy   = tc.true_pos_xy;
                c.tx_power_dBm  = tc.tx_power_dBm;
                c.priority      = priority_map('trusted_fixed');
                c.is_continuing = true;  % trusted 一旦到达就持续
                c.template_idx  = tc.template_idx;
                c.system_knows_position = true;
                c.system_knows_time     = false;
                all_candidates = [all_candidates, c]; %#ok<AGROW>
            end
        end

        % --- 2. 收集 prior 候选 ---
        if ~isempty(PriorCandidatesPerBand(b).candidates)
            for ci = 1:numel(PriorCandidatesPerBand(b).candidates)
                pc = PriorCandidatesPerBand(b).candidates(ci);
                c.instance_id   = pc.instance_id;
                c.template_key  = pc.template_key;
                c.source_type   = pc.source_type;
                c.source_subtype = pc.source_subtype;
                c.true_pos_xy   = pc.true_pos_xy;
                c.tx_power_dBm  = pc.tx_power_dBm;
                c.priority      = priority_map(pc.source_subtype);
                c.is_continuing = pc.is_continuing;
                c.template_idx  = pc.template_idx;
                if strcmp(pc.source_subtype, 'prior_pos_known')
                    c.system_knows_position = true;
                    c.system_knows_time     = true;
                else
                    c.system_knows_position = false;
                    c.system_knows_time     = true;
                end
                all_candidates = [all_candidates, c]; %#ok<AGROW>
            end
        end

        % --- 3. 收集 target 候选 ---
        if ~isempty(TargetCandidatesPerBand(b).candidates)
            for ci = 1:numel(TargetCandidatesPerBand(b).candidates)
                tc = TargetCandidatesPerBand(b).candidates(ci);
                c.instance_id   = tc.instance_id;
                c.template_key  = tc.template_key;
                c.source_type   = tc.source_type;
                c.source_subtype = tc.source_subtype;
                c.true_pos_xy   = tc.true_pos_xy;
                c.tx_power_dBm  = tc.tx_power_dBm;
                c.priority      = priority_map('ordinary_target');
                c.is_continuing = tc.is_continuing;
                c.template_idx  = tc.template_idx;
                c.system_knows_position = false;
                c.system_knows_time     = false;
                all_candidates = [all_candidates, c]; %#ok<AGROW>
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
            % 按复合优先级排序
            % 排序键：[priority DESC, is_continuing DESC, tx_power DESC, template_idx ASC]
            N_cand = numel(all_candidates);
            sort_keys = zeros(N_cand, 4);
            for ci = 1:N_cand
                sort_keys(ci, 1) = all_candidates(ci).priority;
                sort_keys(ci, 2) = double(all_candidates(ci).is_continuing);
                sort_keys(ci, 3) = all_candidates(ci).tx_power_dBm;
                sort_keys(ci, 4) = -all_candidates(ci).template_idx;  % 负号使小编号优先
            end

            % 按列优先级排序（先 priority，再 continuing，再 power，再 idx）
            [~, order] = sortrows(sort_keys, [-1, -2, -3, -4]);
            winner_idx = order(1);
            winner = all_candidates(winner_idx);

            BandWinners(b).has_source = true;
            BandWinners(b).winner     = winner;

            M0State.active_per_band(b).has_source   = true;
            M0State.active_per_band(b).instance_id  = winner.instance_id;
            M0State.active_per_band(b).template_key = winner.template_key;
        end
    end
end
