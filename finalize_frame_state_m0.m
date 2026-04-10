function FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config)
% FINALIZE_FRAME_STATE_M0  组装当前帧的输出结构
%
%   FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config)
%
%   输出 FrameState_t 包含：
%     .frame_id                   - 帧号
%     .time_sec                   - 时间 (s)
%     .active_source_count_unique - 本帧活跃的独立源数量
%     .active_band_count          - 有源的频带数量
%     .active_per_band(b)         - 每频带详细信息

    B = Config.m0.num_bands;
    dt = Config.m0.dt;

    FrameState_t.frame_id = t;
    FrameState_t.time_sec = t * dt;

    active_ids = [];
    active_band_count = 0;

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
            apb.system_knows_position  = w.system_knows_position;
            apb.system_knows_time      = w.system_knows_time;

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
            apb.system_knows_position  = false;
            apb.system_knows_time      = false;

            FrameState_t.active_per_band(b) = apb;
        end
    end

    FrameState_t.active_source_count_unique = numel(active_ids);
    FrameState_t.active_band_count          = active_band_count;
end
