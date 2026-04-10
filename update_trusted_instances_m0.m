function [M0State, TrustedCandidates] = update_trusted_instances_m0(M0State, SourceTemplates, Config, t)
% UPDATE_TRUSTED_INSTANCES_M0  更新宽带 trusted_fixed 源的到达与寿命
%
%   [M0State, TrustedCandidates] = update_trusted_instances_m0(M0State, SourceTemplates, Config, t)
%
%   逻辑：
%     - 对每个 trusted 模板 j：
%       若 inactive → 按泊松到达概率尝试激活，激活时创建新实例
%       若 active   → 寿命减 1，若到 0 则失活
%     - active 时覆盖所有频带 (1:B)
%
%   输出：
%     TrustedCandidates(j).is_active     - 是否活跃
%     TrustedCandidates(j).instance_id   - 实例 ID
%     TrustedCandidates(j).template_key  - 模板 key
%     TrustedCandidates(j).bands_covered - 覆盖频带列表
%     TrustedCandidates(j).true_pos_xy   - 真实位置
%     TrustedCandidates(j).tx_power_dBm  - 发射功率

    dt = Config.m0.dt;
    TrustedNum = Config.m0.trusted.TrustedNum;
    B = Config.m0.num_bands;

    for j = 1:TrustedNum
        tpl = SourceTemplates.trusted(j);

        if ~M0State.trusted_active(j)
            % --- 当前 inactive：尝试到达 ---
            p_arrive = 1 - exp(-tpl.lambda_arrival * dt);
            if rand() < p_arrive
                % 激活：创建新实例
                M0State.trusted_active(j) = true;

                % 采样持续时间
                life = sample_life(tpl.life_mode, tpl.life_param, Config);
                M0State.trusted_life(j) = life;

                % 创建实例
                inst.instance_id          = M0State.next_instance_id;
                inst.template_key         = tpl.template_key;
                inst.source_type          = 'trusted_fixed';
                inst.source_subtype       = 'trusted_fixed';
                inst.band_id_list         = 1:B;
                inst.true_pos_xy          = tpl.fixed_pos_xy;
                inst.tx_power_dBm         = tpl.tx_power_dBm;
                inst.life_remaining       = life;
                inst.is_active            = true;
                inst.system_knows_position = true;
                inst.system_knows_time     = false;

                M0State.next_instance_id = M0State.next_instance_id + 1;

                % 追加到实例池
                if isempty(M0State.instances)
                    M0State.instances = inst;
                else
                    M0State.instances(end+1) = inst;
                end

                M0State.trusted_inst_id(j) = inst.instance_id;
            end
        else
            % --- 当前 active：寿命减 1 ---
            M0State.trusted_life(j) = M0State.trusted_life(j) - 1;

            % 更新实例池中对应实例的 life_remaining
            inst_id = M0State.trusted_inst_id(j);
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == inst_id
                    M0State.instances(k).life_remaining = M0State.trusted_life(j);
                    break;
                end
            end

            if M0State.trusted_life(j) <= 0
                % 失活
                M0State.trusted_active(j) = false;
                for k = 1:numel(M0State.instances)
                    if M0State.instances(k).instance_id == inst_id
                        M0State.instances(k).is_active = false;
                        break;
                    end
                end
            end
        end

        % 构造候选输出
        TrustedCandidates(j).is_active = M0State.trusted_active(j);
        if M0State.trusted_active(j)
            TrustedCandidates(j).instance_id   = M0State.trusted_inst_id(j);
            TrustedCandidates(j).template_key  = tpl.template_key;
            TrustedCandidates(j).bands_covered = 1:B;
            TrustedCandidates(j).true_pos_xy   = tpl.fixed_pos_xy;
            TrustedCandidates(j).tx_power_dBm  = tpl.tx_power_dBm;
            TrustedCandidates(j).template_idx  = j;
        else
            TrustedCandidates(j).instance_id   = 0;
            TrustedCandidates(j).template_key  = '';
            TrustedCandidates(j).bands_covered = [];
            TrustedCandidates(j).true_pos_xy   = [0, 0];
            TrustedCandidates(j).tx_power_dBm  = 0;
            TrustedCandidates(j).template_idx  = j;
        end
    end
end


function life = sample_life(mode, param, Config)
% SAMPLE_LIFE  采样持续帧数
    switch mode
        case 'geom'
            % 几何分布：p_end = param，期望 = 1/param
            life = geornd(param) + 1;  % 至少 1 帧
        case 'exp'
            % 离散指数：mu = param (秒)
            dt = Config.m0.dt;
            life = max(1, ceil(exprnd(param) / dt));
        otherwise
            life = 10;
    end
end
