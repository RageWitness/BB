function [M0State, PriorCandidatesPerBand] = update_prior_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config, t)
% UPDATE_PRIOR_INSTANCES_M0  更新 prior_pos_known 和 prior_time_known 源
%
%   [M0State, PriorCandidatesPerBand] = update_prior_instances_m0(...)
%
%   逻辑：
%     prior_pos_known  — 按 schedule/periodic 判断当前帧是否应出现
%     prior_time_known — 按 schedule/periodic 判断；每次激活起点采样真实位置
%
%   输出：
%     PriorCandidatesPerBand(b).candidates — 该频带所有活跃 prior 候选

    B = Config.m0.num_bands;

    % 初始化输出
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
                % 查找是否已有活跃实例
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);

                if inst_idx == 0
                    % 创建新实例
                    inst.instance_id          = M0State.next_instance_id;
                    inst.template_key         = tpl.template_key;
                    inst.source_type          = 'prior';
                    inst.source_subtype       = 'prior_pos_known';
                    inst.band_id_list         = b;
                    inst.true_pos_xy          = tpl.fixed_pos_xy;
                    inst.tx_power_dBm         = tpl.tx_power_dBm;
                    inst.life_remaining       = -1;  % schedule 控制，不用 life
                    inst.is_active            = true;
                    inst.system_knows_position = true;
                    inst.system_knows_time     = true;

                    M0State.next_instance_id = M0State.next_instance_id + 1;
                    if isempty(M0State.instances)
                        M0State.instances = inst;
                    else
                        M0State.instances(end+1) = inst;
                    end
                    inst_idx = numel(M0State.instances);
                else
                    M0State.instances(inst_idx).is_active = true;
                end

                % 添加为候选
                cand.instance_id   = M0State.instances(inst_idx).instance_id;
                cand.template_key  = tpl.template_key;
                cand.source_type   = 'prior';
                cand.source_subtype = 'prior_pos_known';
                cand.true_pos_xy   = tpl.fixed_pos_xy;
                cand.tx_power_dBm  = tpl.tx_power_dBm;
                cand.template_idx  = s;
                cand.is_continuing = (M0State.instances(inst_idx).life_remaining == -1);

                PriorCandidatesPerBand(b).candidates = [ ...
                    PriorCandidatesPerBand(b).candidates, cand];
            else
                % 如有实例则标记失活
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
                    % 新激活起点 —— 采样真实位置
                    pos_idx = randi(GridValid.Nvalid);
                    pos_xy  = GridValid.xy(pos_idx, :);

                    inst.instance_id          = M0State.next_instance_id;
                    inst.template_key         = tpl.template_key;
                    inst.source_type          = 'prior';
                    inst.source_subtype       = 'prior_time_known';
                    inst.band_id_list         = b;
                    inst.true_pos_xy          = pos_xy;
                    inst.tx_power_dBm         = tpl.tx_power_dBm;
                    inst.life_remaining       = -1;
                    inst.is_active            = true;
                    inst.system_knows_position = false;
                    inst.system_knows_time     = true;

                    M0State.next_instance_id = M0State.next_instance_id + 1;
                    if isempty(M0State.instances)
                        M0State.instances = inst;
                    else
                        M0State.instances(end+1) = inst;
                    end
                    inst_idx = numel(M0State.instances);
                else
                    % 已在激活期，沿用本次事件位置
                    M0State.instances(inst_idx).is_active = true;
                end

                cand.instance_id   = M0State.instances(inst_idx).instance_id;
                cand.template_key  = tpl.template_key;
                cand.source_type   = 'prior';
                cand.source_subtype = 'prior_time_known';
                cand.true_pos_xy   = M0State.instances(inst_idx).true_pos_xy;
                cand.tx_power_dBm  = tpl.tx_power_dBm;
                cand.template_idx  = q;
                cand.is_continuing = true;

                PriorCandidatesPerBand(b).candidates = [ ...
                    PriorCandidatesPerBand(b).candidates, cand];
            else
                % schedule 关闭时，删除该实例以便下次重新采样位置
                inst_idx = find_instance_by_key(M0State.instances, tpl.template_key);
                if inst_idx > 0
                    M0State.instances(inst_idx).is_active = false;
                    % 清除实例，使下次激活时重新采样位置
                    M0State.instances(inst_idx).life_remaining = 0;
                end
            end
        end
    end
end


function is_on = check_schedule(tpl, t)
% CHECK_SCHEDULE  检查模板在帧 t 是否应活跃
    if strcmp(tpl.schedule_mode, 'periodic')
        phase = mod(t - 1 - tpl.phase_frames, tpl.period_frames);
        is_on = phase < tpl.duration_frames;
    else
        % timetable 模式
        is_on = ismember(t, tpl.schedule_set);
    end
end


function idx = find_instance_by_key(instances, key)
% FIND_INSTANCE_BY_KEY  在实例池中按 template_key 查找活跃实例
%   返回 0 表示未找到
    idx = 0;
    if isempty(instances)
        return;
    end
    for k = 1:numel(instances)
        if strcmp(instances(k).template_key, key) && instances(k).is_active
            idx = k;
            return;
        end
    end
    % 如果没找到活跃的，查找 life_remaining == -1 的（prior schedule 控制的）
    for k = 1:numel(instances)
        if strcmp(instances(k).template_key, key) && instances(k).life_remaining == -1
            idx = k;
            return;
        end
    end
end
