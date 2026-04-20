function [M0State, TargetCandidatesPerBand] = update_target_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config)
% UPDATE_TARGET_INSTANCES_M0  待定位源的到达、位置、功率、寿命（每频带一路）
%
%   输出 TargetCandidatesPerBand(b).candidates 为 0 或 1 个候选（新 SourceEvent 片段）

    B = Config.m0.num_bands;
    dt = Config.m0.dt;

    for b = 1:B
        TargetCandidatesPerBand(b).candidates = [];
    end

    for b = 1:B
        was_active = M0State.target_active_band(b);

        if was_active
            M0State.target_life_band(b) = M0State.target_life_band(b) - 1;
            if M0State.target_life_band(b) <= 0
                M0State.target_active_band(b) = false;
                M0State.target_inst_band(b) = 0;
            end
        end

        if M0State.target_active_band(b)
            % 从实例池拿出对应实例
            inst_id = M0State.target_inst_band(b);
            inst = find_instance(M0State, inst_id);
            if ~isempty(inst)
                TargetCandidatesPerBand(b).candidates = build_target_candidate( ...
                    inst, b, 0, true);
            end
        else
            % 尝试到达（从模板池中选一个候选）
            Nocoop_b = Config.m0.target.Nocoop(b);
            for n = 1:Nocoop_b
                tpl = SourceTemplates.target(b, n);
                p = 1 - exp(-tpl.lambda_arrival * dt);
                if rand() < p
                    pos_xy = sample_position(GridValid, tpl.position_mode);
                    power  = sample_power(tpl.tx_power_range_dBm);
                    life   = sample_target_life(tpl.life_mode, tpl.life_param, Config);

                    inst = struct();
                    inst.instance_id      = M0State.next_instance_id;
                    inst.template_key     = tpl.template_key;
                    inst.source_type      = 'target';
                    inst.band_list        = b;
                    inst.true_position    = pos_xy;
                    inst.true_tx_power    = power;
                    inst.tx_power_by_band = [];
                    inst.location_prior   = struct('type','none','value', []);
                    inst.power_prior      = struct('type','none','value', []);
                    inst.life_remaining   = life;
                    inst.is_active        = true;

                    M0State.next_instance_id = M0State.next_instance_id + 1;
                    M0State.instances = append_instance(M0State.instances, inst);

                    M0State.target_active_band(b) = true;
                    M0State.target_life_band(b)   = life;
                    M0State.target_inst_band(b)   = inst.instance_id;

                    TargetCandidatesPerBand(b).candidates = build_target_candidate( ...
                        inst, b, n, false);
                    break;
                end
            end
        end
    end
end


function cand = build_target_candidate(inst, b, template_idx, is_continuing)
    cand.instance_id      = inst.instance_id;
    cand.template_key     = inst.template_key;
    cand.source_type      = 'target';
    cand.source_subtype   = 'target';
    cand.band_list        = b;
    cand.band_id          = b;
    cand.true_position    = inst.true_position;
    cand.true_tx_power    = inst.true_tx_power;
    cand.tx_power_by_band = [];
    cand.is_continuing    = is_continuing;
    cand.template_idx     = template_idx;
    cand.is_persistent    = false;
    cand.location_prior   = inst.location_prior;
    cand.power_prior      = inst.power_prior;
end


function inst = find_instance(M0State, inst_id)
    inst = [];
    for k = 1:numel(M0State.instances)
        if M0State.instances(k).instance_id == inst_id
            inst = M0State.instances(k);
            return;
        end
    end
end


function arr = append_instance(arr, inst)
    if isempty(arr)
        arr = inst;
    else
        arr(end+1) = inst;
    end
end


function pos_xy = sample_position(GridValid, mode)
    switch mode
        case 'uniform'
            idx = randi(GridValid.Nvalid);
            pos_xy = GridValid.xy(idx, :);
        case 'hotspot'
            centers = [75, 75; 150, 225; 225, 150];
            sigma = 40;
            c_idx = randi(size(centers, 1));
            raw = centers(c_idx, :) + sigma * randn(1, 2);
            dists = sum((GridValid.xy - raw).^2, 2);
            [~, idx] = min(dists);
            pos_xy = GridValid.xy(idx, :);
        otherwise
            idx = randi(GridValid.Nvalid);
            pos_xy = GridValid.xy(idx, :);
    end
end


function power = sample_power(range_dBm)
    power = range_dBm(1) + (range_dBm(2) - range_dBm(1)) * rand();
end


function life = sample_target_life(mode, param, Config)
    switch mode
        case 'geom'
            life = geornd(param) + 1;
        case 'uniform'
            L_range = Config.m0.target.life_range;
            life = randi([L_range(1), L_range(2)]);
        case 'exp'
            life = max(1, ceil(exprnd(param) / Config.m0.dt));
        otherwise
            life = 10;
    end
end
