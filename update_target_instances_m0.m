function [M0State, TargetCandidatesPerBand] = update_target_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config, t)
% UPDATE_TARGET_INSTANCES_M0  更新 ordinary_target 的到达、位置、功率、寿命
%
%   [M0State, TargetCandidatesPerBand] = update_target_instances_m0(...)
%
%   逻辑（对每个频带 b）：
%     若当前没有持续中的 ordinary target → 按泊松概率尝试产生新目标
%       新目标：采样位置、功率、持续时间
%     若有持续中的 ordinary target → 寿命减 1，到 0 则失活
%
%   输出：
%     TargetCandidatesPerBand(b).candidates — 该频带活跃的 ordinary target 候选

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

            % 更新实例池
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
            % 仍然活跃，构造候选
            inst_id = M0State.target_inst_band(b);
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == inst_id
                    cand.instance_id    = inst_id;
                    cand.template_key   = M0State.instances(k).template_key;
                    cand.source_type    = 'ordinary_target';
                    cand.source_subtype = 'ordinary_target';
                    cand.true_pos_xy    = M0State.instances(k).true_pos_xy;
                    cand.tx_power_dBm   = M0State.instances(k).tx_power_dBm;
                    cand.template_idx   = 0;
                    cand.is_continuing  = true;

                    TargetCandidatesPerBand(b).candidates = cand;
                    break;
                end
            end
        else
            % --- 当前无目标：尝试泊松到达 ---
            % 遍历该频带的所有模板，选第一个成功到达的
            Nocoop_b = Config.m0.target.Nocoop(b);
            for n = 1:Nocoop_b
                tpl = SourceTemplates.target(b, n);
                p_arrive = 1 - exp(-tpl.lambda_arrival * dt);

                if rand() < p_arrive
                    % 到达！采样位置、功率、持续时间
                    pos_xy = sample_position(GridValid, tpl.position_mode);
                    power  = sample_power(tpl.tx_power_range_dBm);
                    life   = sample_target_life(tpl.life_mode, tpl.life_param, Config);

                    % 创建实例
                    inst.instance_id          = M0State.next_instance_id;
                    inst.template_key         = tpl.template_key;
                    inst.source_type          = 'ordinary_target';
                    inst.source_subtype       = 'ordinary_target';
                    inst.band_id_list         = b;
                    inst.true_pos_xy          = pos_xy;
                    inst.tx_power_dBm         = power;
                    inst.life_remaining       = life;
                    inst.is_active            = true;
                    inst.system_knows_position = false;
                    inst.system_knows_time     = false;

                    M0State.next_instance_id = M0State.next_instance_id + 1;
                    if isempty(M0State.instances)
                        M0State.instances = inst;
                    else
                        M0State.instances(end+1) = inst;
                    end

                    % 记录频带追踪信息
                    M0State.target_active_band(b) = true;
                    M0State.target_life_band(b)   = life;
                    M0State.target_inst_band(b)   = inst.instance_id;

                    % 构造候选
                    cand.instance_id    = inst.instance_id;
                    cand.template_key   = tpl.template_key;
                    cand.source_type    = 'ordinary_target';
                    cand.source_subtype = 'ordinary_target';
                    cand.true_pos_xy    = pos_xy;
                    cand.tx_power_dBm   = power;
                    cand.template_idx   = n;
                    cand.is_continuing  = false;

                    TargetCandidatesPerBand(b).candidates = cand;

                    break;  % 该频带只产生 1 个新目标
                end
            end
        end
    end
end


function pos_xy = sample_position(GridValid, mode)
% SAMPLE_POSITION  从有效网格采样位置
    switch mode
        case 'uniform'
            idx = randi(GridValid.Nvalid);
            pos_xy = GridValid.xy(idx, :);
        case 'hotspot'
            % 热点模式：以几个中心点为高斯混合采样，然后 snap 到最近网格点
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
% SAMPLE_POWER  从功率范围均匀采样
    power = range_dBm(1) + (range_dBm(2) - range_dBm(1)) * rand();
end


function life = sample_target_life(mode, param, Config)
% SAMPLE_TARGET_LIFE  采样 ordinary target 持续帧数
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
