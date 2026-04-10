function FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config)
% FINALIZE_FRAME_STATE_M0  组装当前帧的输出结构（含先验对象）
%
%   FrameState_t.active_per_band(b) 现在除了真值信息外，还携带：
%     position_prior, time_prior, power_nominal_dBm, power_range_dBm,
%     power_stability_level, credibility_prior, upgrade_potential

    B = Config.m0.num_bands;
    dt = Config.m0.dt;

    FrameState_t.frame_id = t;
    FrameState_t.time_sec = t * dt;

    active_ids = [];
    active_band_count = 0;

    % 空先验对象（用于无源频带）
    empty_pos_prior  = struct('mode','','candidate_positions',[],'match_radius',[]);
    empty_time_prior = struct('mode','','schedule',[],'period_frames',[], ...
                              'duration_frames',[],'phase_frames',[]);

    for b = 1:B
        if BandWinners(b).has_source
            w = BandWinners(b).winner;

            apb.has_source             = true;
            apb.source_type            = w.source_type;
            apb.source_subtype         = w.source_subtype;
            apb.template_key           = w.template_key;
            apb.instance_id            = w.instance_id;
            apb.band_id                = b;
            apb.true_pos_xy            = w.true_pos_xy;
            apb.tx_power_dBm           = w.tx_power_dBm;

            % 查找 life_remaining
            apb.life_remaining = 0;
            for k = 1:numel(M0State.instances)
                if M0State.instances(k).instance_id == w.instance_id
                    apb.life_remaining = M0State.instances(k).life_remaining;
                    break;
                end
            end

            % 标记
            apb.is_calibrator          = strcmp(w.source_type, 'trusted_fixed') || ...
                                         strcmp(w.source_subtype, 'prior_pos_known');
            apb.is_localization_target = strcmp(w.source_type, 'ordinary_target');

            % --- 先验对象（从 winner 直接传递） ---
            apb.position_prior        = w.position_prior;
            apb.time_prior            = w.time_prior;
            apb.power_nominal_dBm     = w.power_nominal_dBm;
            apb.power_range_dBm       = w.power_range_dBm;
            apb.power_stability_level = w.power_stability_level;
            apb.credibility_prior     = w.credibility_prior;
            apb.upgrade_potential     = w.upgrade_potential;

            FrameState_t.active_per_band(b) = apb;

            active_band_count = active_band_count + 1;
            if ~ismember(w.instance_id, active_ids)
                active_ids(end+1) = w.instance_id; %#ok<AGROW>
            end
        else
            apb.has_source             = false;
            apb.source_type            = '';
            apb.source_subtype         = '';
            apb.template_key           = '';
            apb.instance_id            = 0;
            apb.band_id                = b;
            apb.true_pos_xy            = [0, 0];
            apb.tx_power_dBm           = 0;
            apb.life_remaining         = 0;
            apb.is_calibrator          = false;
            apb.is_localization_target = false;
            apb.position_prior         = empty_pos_prior;
            apb.time_prior             = empty_time_prior;
            apb.power_nominal_dBm     = 0;
            apb.power_range_dBm       = [0, 0];
            apb.power_stability_level = 0;
            apb.credibility_prior     = 0;
            apb.upgrade_potential     = 'none';

            FrameState_t.active_per_band(b) = apb;
        end
    end

    FrameState_t.active_source_count_unique = numel(active_ids);
    FrameState_t.active_band_count          = active_band_count;
end
