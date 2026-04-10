function [M0State, TargetCandidatesPerBand] = update_target_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config, t)
% UPDATE_TARGET_INSTANCES_M0  更新 ordinary_target 的到达、位置、功率、寿命
%
%   实例创建时绑定模板的 position_prior, time_prior, upgrade_potential 等。
%   候选输出同样携带完整先验对象。

    B = Config.m0.num_bands;
    dt = Config.m0.dt;

    for b = 1:B
        TargetCandidatesPerBand(b).candidates = [];
    end

    for b = 1:B
        if M0State.target_active_band(b)
            % --- 当前有持续中的目标：寿命减 1 ---
            M0State.target_life_band(b) = M0State.target_life_band(b) - 1;

            inst_id = M0State.target_inst_band(b);
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == inst_id
                    M0State.instances(k).life_remaining = M0State.target_life_band(b);
                    if M0State.target_life_band(b) <= 0
                        M0State.instances(k).is_active = false;
                        M0State.target_active_band(b) = false;
                        M0State.target_inst_band(b) = 0;
                    end
                    break;
                end
            end
        end

        if M0State.target_active_band(b)
            % 仍然活跃，构造候选（从实例读取全部先验字段）
            inst_id = M0State.target_inst_band(b);
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == inst_id
                    cand = build_target_candidate(M0State.instances(k), 0, true);
                    TargetCandidatesPerBand(b).candidates = cand;
                    break;
                end
            end
        else
            % --- 当前无目标：尝试泊松到达 ---
            Nocoop_b = Config.m0.target.Nocoop(b);
            for n = 1:Nocoop_b
                tpl = SourceTemplates.target(b, n);
                p_arrive = 1 - exp(-tpl.lambda_arrival * dt);

                if rand() < p_arrive
                    pos_xy = sample_position(GridValid, tpl.position_mode);
                    power  = sample_power(tpl.tx_power_range_dBm);
                    life   = sample_target_life(tpl.life_mode, tpl.life_param, Config);

                    % 创建实例（含完整先验对象）
                    inst.instance_id           = M0State.next_instance_id;
                    inst.template_key          = tpl.template_key;
                    inst.source_type           = 'ordinary_target';
                    inst.source_subtype        = 'ordinary_target';
                    inst.band_id_list          = b;
                    inst.true_pos_xy           = pos_xy;
                    inst.tx_power_dBm          = power;
                    inst.life_remaining        = life;
                    inst.is_active             = true;
                    inst.position_prior        = tpl.position_prior;   % mode='unknown'
                    inst.time_prior            = tpl.time_prior;       % mode='unknown'
                    inst.power_nominal_dBm     = tpl.power_nominal_dBm;
                    inst.power_range_dBm       = tpl.power_range_dBm;
                    inst.power_stability_level = tpl.power_stability_level;
                    inst.credibility_prior     = tpl.credibility_prior;
                    inst.upgrade_potential     = tpl.upgrade_potential;

                    M0State.next_instance_id = M0State.next_instance_id + 1;
                    if isempty(M0State.instances)
                        M0State.instances = inst;
                    else
                        M0State.instances(end+1) = inst; %#ok<AGROW>
                    end

                    M0State.target_active_band(b) = true;
                    M0State.target_life_band(b)   = life;
                    M0State.target_inst_band(b)   = inst.instance_id;

                    cand = build_target_candidate(inst, n, false);
                    TargetCandidatesPerBand(b).candidates = cand;

                    break;
                end
            end
        end
    end
end


function cand = build_target_candidate(inst, template_idx, is_continuing)
% BUILD_TARGET_CANDIDATE  从实例构造带先验对象的 target 候选
    cand.instance_id           = inst.instance_id;
    cand.template_key          = inst.template_key;
    cand.source_type           = 'ordinary_target';
    cand.source_subtype        = 'ordinary_target';
    cand.true_pos_xy           = inst.true_pos_xy;
    cand.tx_power_dBm          = inst.tx_power_dBm;
    cand.template_idx          = template_idx;
    cand.is_continuing         = is_continuing;
    cand.position_prior        = inst.position_prior;
    cand.time_prior            = inst.time_prior;
    cand.power_nominal_dBm     = inst.power_nominal_dBm;
    cand.power_range_dBm       = inst.power_range_dBm;
    cand.power_stability_level = inst.power_stability_level;
    cand.credibility_prior     = inst.credibility_prior;
    cand.upgrade_potential     = inst.upgrade_potential;
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
            dt = Config.m0.dt;
            life = max(1, ceil(exprnd(param) / dt));
        otherwise
            life = 10;
    end
end
