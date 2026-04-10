function [M0State, TrustedCandidates] = update_trusted_instances_m0(M0State, SourceTemplates, Config, t)
% UPDATE_TRUSTED_INSTANCES_M0  更新宽带 trusted_fixed 源的到达与寿命
%
%   创建实例时绑定模板的 position_prior, time_prior 等先验对象。
%   候选输出同样携带先验对象，供裁决与帧输出使用。

    dt = Config.m0.dt;
    TrustedNum = Config.m0.trusted.TrustedNum;
    B = Config.m0.num_bands;

    for j = 1:TrustedNum
        tpl = SourceTemplates.trusted(j);

        if ~M0State.trusted_active(j)
            % --- 当前 inactive：尝试到达 ---
            p_arrive = 1 - exp(-tpl.lambda_arrival * dt);
            if rand() < p_arrive
                M0State.trusted_active(j) = true;

                life = sample_life(tpl.life_mode, tpl.life_param, Config);
                M0State.trusted_life(j) = life;

                % 创建实例（含先验对象）
                inst.instance_id           = M0State.next_instance_id;
                inst.template_key          = tpl.template_key;
                inst.source_type           = 'trusted_fixed';
                inst.source_subtype        = 'trusted_fixed';
                inst.band_id_list          = 1:B;
                inst.true_pos_xy           = tpl.fixed_pos_xy;
                inst.tx_power_dBm          = tpl.tx_power_dBm;
                inst.life_remaining        = life;
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
                    M0State.instances(end+1) = inst;
                end

                M0State.trusted_inst_id(j) = inst.instance_id;
            end
        else
            % --- 当前 active：寿命减 1 ---
            M0State.trusted_life(j) = M0State.trusted_life(j) - 1;

            inst_id = M0State.trusted_inst_id(j);
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == inst_id
                    M0State.instances(k).life_remaining = M0State.trusted_life(j);
                    break;
                end
            end

            if M0State.trusted_life(j) <= 0
                M0State.trusted_active(j) = false;
                for k = 1:numel(M0State.instances)
                    if M0State.instances(k).instance_id == inst_id
                        M0State.instances(k).is_active = false;
                        break;
                    end
                end
            end
        end

        % 构造候选输出（含先验对象）
        TrustedCandidates(j).is_active = M0State.trusted_active(j);
        if M0State.trusted_active(j)
            TrustedCandidates(j).instance_id           = M0State.trusted_inst_id(j);
            TrustedCandidates(j).template_key          = tpl.template_key;
            TrustedCandidates(j).bands_covered         = 1:B;
            TrustedCandidates(j).true_pos_xy           = tpl.fixed_pos_xy;
            TrustedCandidates(j).tx_power_dBm          = tpl.tx_power_dBm;
            TrustedCandidates(j).template_idx          = j;
            TrustedCandidates(j).position_prior        = tpl.position_prior;
            TrustedCandidates(j).time_prior            = tpl.time_prior;
            TrustedCandidates(j).power_nominal_dBm     = tpl.power_nominal_dBm;
            TrustedCandidates(j).power_range_dBm       = tpl.power_range_dBm;
            TrustedCandidates(j).power_stability_level = tpl.power_stability_level;
            TrustedCandidates(j).credibility_prior     = tpl.credibility_prior;
            TrustedCandidates(j).upgrade_potential     = tpl.upgrade_potential;
        else
            TrustedCandidates(j).instance_id           = 0;
            TrustedCandidates(j).template_key          = '';
            TrustedCandidates(j).bands_covered         = [];
            TrustedCandidates(j).true_pos_xy           = [0, 0];
            TrustedCandidates(j).tx_power_dBm          = 0;
            TrustedCandidates(j).template_idx          = j;
            TrustedCandidates(j).position_prior        = struct('mode','','candidate_positions',[],'match_radius',[]);
            TrustedCandidates(j).time_prior            = struct('mode','','schedule',[],'period_frames',[],'duration_frames',[],'phase_frames',[]);
            TrustedCandidates(j).power_nominal_dBm     = 0;
            TrustedCandidates(j).power_range_dBm       = [0, 0];
            TrustedCandidates(j).power_stability_level = 0;
            TrustedCandidates(j).credibility_prior     = 0;
            TrustedCandidates(j).upgrade_potential     = 'none';
        end
    end
end


function life = sample_life(mode, param, Config)
% SAMPLE_LIFE  采样持续帧数
    switch mode
        case 'geom'
            life = geornd(param) + 1;
        case 'exp'
            dt = Config.m0.dt;
            life = max(1, ceil(exprnd(param) / dt));
        otherwise
            life = 10;
    end
end
