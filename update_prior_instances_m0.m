function [M0State, PriorCandidatesPerBand] = update_prior_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config, t)
% UPDATE_PRIOR_INSTANCES_M0  更新 prior_pos_known 和 prior_time_known 源
%
%   实例创建时绑定 position_prior / time_prior 等先验对象。
%   候选输出同样携带先验对象。

    B = Config.m0.num_bands;

    for b = 1:B
        PriorCandidatesPerBand(b).candidates = [];
    end

    %% ===== prior_pos_known =====
    for b = 1:B
        Sprior_b = Config.m0.prior.Sprior(b);
        for s = 1:Sprior_b
            tpl = SourceTemplates.prior_pos(b, s);
            is_on = check_schedule(tpl, t);

            if is_on
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);

                if inst_idx == 0
                    inst.instance_id           = M0State.next_instance_id;
                    inst.template_key          = tpl.template_key;
                    inst.source_type           = 'prior';
                    inst.source_subtype        = 'prior_pos_known';
                    inst.band_id_list          = b;
                    inst.true_pos_xy           = tpl.fixed_pos_xy;
                    inst.tx_power_dBm          = tpl.tx_power_dBm;
                    inst.life_remaining        = -1;
                    inst.is_active             = true;
                    inst.position_prior        = tpl.position_prior;
                    inst.time_prior            = tpl.time_prior;
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
                    inst_idx = numel(M0State.instances);
                else
                    M0State.instances(inst_idx).is_active = true;
                end

                % 构造候选
                cand = build_prior_candidate(M0State.instances(inst_idx), tpl, s);
                PriorCandidatesPerBand(b).candidates = [ ...
                    PriorCandidatesPerBand(b).candidates, cand];
            else
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);
                if inst_idx > 0
                    M0State.instances(inst_idx).is_active = false;
                end
            end
        end
    end

    %% ===== prior_time_known =====
    for b = 1:B
        Tprior_b = Config.m0.prior.Tprior(b);
        for q = 1:Tprior_b
            tpl = SourceTemplates.prior_time(b, q);
            is_on = check_schedule(tpl, t);

            if is_on
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);

                if inst_idx == 0
                    % 新激活起点 —— 采样真实位置（系统不知道）
                    pos_idx = randi(GridValid.Nvalid);
                    pos_xy  = GridValid.xy(pos_idx, :);

                    inst.instance_id           = M0State.next_instance_id;
                    inst.template_key          = tpl.template_key;
                    inst.source_type           = 'prior';
                    inst.source_subtype        = 'prior_time_known';
                    inst.band_id_list          = b;
                    inst.true_pos_xy           = pos_xy;
                    inst.tx_power_dBm          = tpl.tx_power_dBm;
                    inst.life_remaining        = -1;
                    inst.is_active             = true;
                    % 位置先验：系统只知道时间，不知道位置
                    inst.position_prior        = tpl.position_prior;  % mode='time_known_position_unknown'
                    inst.time_prior            = tpl.time_prior;
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
                    inst_idx = numel(M0State.instances);
                else
                    M0State.instances(inst_idx).is_active = true;
                end

                cand = build_prior_candidate(M0State.instances(inst_idx), tpl, q);
                PriorCandidatesPerBand(b).candidates = [ ...
                    PriorCandidatesPerBand(b).candidates, cand];
            else
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);
                if inst_idx > 0
                    M0State.instances(inst_idx).is_active = false;
                    M0State.instances(inst_idx).life_remaining = 0;
                end
            end
        end
    end
end


function cand = build_prior_candidate(inst, tpl, template_idx)
% BUILD_PRIOR_CANDIDATE  从实例和模板构造带先验对象的候选
    cand.instance_id           = inst.instance_id;
    cand.template_key          = inst.template_key;
    cand.source_type           = inst.source_type;
    cand.source_subtype        = inst.source_subtype;
    cand.true_pos_xy           = inst.true_pos_xy;
    cand.tx_power_dBm          = inst.tx_power_dBm;
    cand.template_idx          = template_idx;
    cand.is_continuing         = (inst.life_remaining == -1);
    cand.position_prior        = inst.position_prior;
    cand.time_prior            = inst.time_prior;
    cand.power_nominal_dBm     = inst.power_nominal_dBm;
    cand.power_range_dBm       = inst.power_range_dBm;
    cand.power_stability_level = inst.power_stability_level;
    cand.credibility_prior     = inst.credibility_prior;
    cand.upgrade_potential     = inst.upgrade_potential;
end


function is_on = check_schedule(tpl, t)
% CHECK_SCHEDULE  检查模板在帧 t 是否应活跃
    if strcmp(tpl.schedule_mode, 'periodic')
        phase = mod(t - 1 - tpl.phase_frames, tpl.period_frames);
        is_on = phase < tpl.duration_frames;
    else
        is_on = ismember(t, tpl.schedule_set);
    end
end


function idx = find_instance_by_key(instances, key)
% FIND_INSTANCE_BY_KEY  按 template_key 查找实例
    idx = 0;
    if isempty(instances), return; end
    for k = 1:numel(instances)
        if strcmp(instances(k).template_key, key) && instances(k).is_active
            idx = k; return;
        end
    end
    for k = 1:numel(instances)
        if strcmp(instances(k).template_key, key) && instances(k).life_remaining == -1
            idx = k; return;
        end
    end
end
