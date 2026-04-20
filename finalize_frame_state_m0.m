function FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config)
% FINALIZE_FRAME_STATE_M0  组装当前帧的统一 SourceEvent 输出
%
%   FrameState_t.active_per_band(b) 字段（统一 SourceEvent 结构）：
%     source_type, priority, band_id, band_list,
%     true_position, true_tx_power, tx_power_by_band,
%     location_prior(.type/.value), power_prior(.type/.value),
%     timestamp, is_persistent, is_selected_final
%   兼容字段（过渡期保留给 M1 / 现有下游代码）：
%     has_source, true_pos_xy, tx_power_dBm,
%     is_calibrator, is_localization_target, template_key, instance_id, life_remaining

    B  = Config.m0.num_bands;
    dt = Config.m0.dt;

    FrameState_t.frame_id = t;
    FrameState_t.time_sec = t * dt;

    active_ids = [];
    active_band_count = 0;

    for b = 1:B
        if BandWinners(b).has_source
            w = BandWinners(b).winner;

            % --- 新 SourceEvent 字段 ---
            apb.source_type        = w.source_type;
            apb.priority           = w.priority;
            apb.band_id            = b;
            if isempty(w.band_list)
                apb.band_list = b;
            else
                apb.band_list = w.band_list;
            end
            apb.true_position      = w.true_position;
            apb.true_tx_power      = w.true_tx_power;
            apb.tx_power_by_band   = w.tx_power_by_band;
            apb.location_prior     = w.location_prior;
            apb.power_prior        = w.power_prior;
            apb.timestamp          = t;
            apb.is_persistent      = w.is_persistent;
            apb.is_selected_final  = true;

            % --- 兼容字段 ---
            apb.has_source             = true;
            apb.source_subtype         = w.source_type;
            apb.template_key           = w.template_key;
            apb.instance_id            = w.instance_id;
            apb.true_pos_xy            = w.true_position;
            apb.tx_power_dBm           = w.true_tx_power;
            apb.is_calibrator          = strcmp(w.source_type, 'broadband_cal') || ...
                                         strcmp(w.source_type, 'persistent_cal');
            apb.is_localization_target = strcmp(w.source_type, 'target');

            % life_remaining（查找实例池；持续源固定 -1）
            apb.life_remaining = -1;
            if w.instance_id > 0
                for k = 1:numel(M0State.instances)
                    if M0State.instances(k).instance_id == w.instance_id
                        apb.life_remaining = M0State.instances(k).life_remaining;
                        break;
                    end
                end
            end

            FrameState_t.active_per_band(b) = apb;

            active_band_count = active_band_count + 1;
            if w.instance_id > 0 && ~ismember(w.instance_id, active_ids)
                active_ids(end+1) = w.instance_id; %#ok<AGROW>
            end
        else
            FrameState_t.active_per_band(b) = empty_event(b, t);
        end
    end

    FrameState_t.active_source_count_unique = numel(active_ids);
    FrameState_t.active_band_count          = active_band_count;
end


function apb = empty_event(b, t)
    apb.source_type        = '';
    apb.priority           = 0;
    apb.band_id            = b;
    apb.band_list          = b;
    apb.true_position      = [0, 0];
    apb.true_tx_power      = 0;
    apb.tx_power_by_band   = [];
    apb.location_prior     = struct('type','none','value',[]);
    apb.power_prior        = struct('type','none','value',[]);
    apb.timestamp          = t;
    apb.is_persistent      = false;
    apb.is_selected_final  = false;

    apb.has_source             = false;
    apb.source_subtype         = '';
    apb.template_key           = '';
    apb.instance_id            = 0;
    apb.true_pos_xy            = [0, 0];
    apb.tx_power_dBm           = 0;
    apb.is_calibrator          = false;
    apb.is_localization_target = false;
    apb.life_remaining         = 0;
end
