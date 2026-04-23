function [M0State, Candidates] = update_broadband_cal_instances_m0( ...
    M0State, SourceTemplates, Config)
% UPDATE_BROADBAND_CAL_INSTANCES_M0  宽带标校源到达与寿命
%
%   支持两种模式：
%     'poisson' — Poisson 到达 + geometric/exp 寿命（原逻辑）
%     'manual'  — 用户指定帧区间，区间内 active，区间外 inactive

    dt = Config.m0.dt;
    t  = M0State.frame_id;
    BC_count = numel(SourceTemplates.broadband_cal);

    Candidates = [];

    for j = 1:BC_count
        tpl = SourceTemplates.broadband_cal(j);
        was_active = M0State.broadband_active(j);

        is_manual = isfield(tpl, 'schedule_mode') && strcmp(tpl.schedule_mode, 'manual');

        if is_manual
            % --- 手动调度：帧区间判定 ---
            in_range = (t >= tpl.frame_range(1)) && (t <= tpl.frame_range(2));

            if in_range && ~was_active
                M0State.broadband_active(j) = true;
                M0State.broadband_life(j)   = tpl.frame_range(2) - t + 1;
                M0State.broadband_inst_id(j) = M0State.next_instance_id;
                M0State.next_instance_id = M0State.next_instance_id + 1;
            elseif in_range && was_active
                M0State.broadband_life(j) = tpl.frame_range(2) - t + 1;
            elseif ~in_range && was_active
                M0State.broadband_active(j) = false;
                M0State.broadband_inst_id(j) = 0;
            end
        else
            % --- Poisson 到达 + 随机寿命（原逻辑）---
            if ~was_active
                p = 1 - exp(-tpl.lambda_arrival * dt);
                if rand() < p
                    M0State.broadband_active(j) = true;
                    M0State.broadband_life(j) = sample_life(tpl.life_mode, tpl.life_param, Config);
                    M0State.broadband_inst_id(j) = M0State.next_instance_id;
                    M0State.next_instance_id = M0State.next_instance_id + 1;
                end
            else
                M0State.broadband_life(j) = M0State.broadband_life(j) - 1;
                if M0State.broadband_life(j) <= 0
                    M0State.broadband_active(j) = false;
                    M0State.broadband_inst_id(j) = 0;
                end
            end
        end

        if M0State.broadband_active(j)
            c = struct();
            c.instance_id      = M0State.broadband_inst_id(j);
            c.template_key     = tpl.template_key;
            c.source_type      = 'broadband_cal';
            c.source_subtype   = 'broadband_cal';
            c.band_list        = tpl.band_list;
            c.band_id          = 0;
            c.true_position    = tpl.fixed_pos_xy;
            c.true_tx_power    = tpl.tx_power_dBm;
            c.tx_power_by_band = tpl.tx_power_by_band;
            c.is_continuing    = was_active;
            c.template_idx     = j;
            c.is_persistent    = false;
            c.location_prior   = struct('type','exact','value', tpl.fixed_pos_xy);
            c.power_prior      = struct('type','exact','value', tpl.tx_power_dBm);

            if isempty(Candidates)
                Candidates = c;
            else
                Candidates(end+1) = c; %#ok<AGROW>
            end
        end
    end
end


function life = sample_life(mode, param, Config)
    switch mode
        case 'geom'
            life = geornd(param) + 1;
        case 'exp'
            life = max(1, ceil(exprnd(param) / Config.m0.dt));
        otherwise
            life = 10;
    end
end
