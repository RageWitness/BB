function M0Logs = append_truth_logs_m0(M0Logs, FrameState_t, Config)
% APPEND_TRUTH_LOGS_M0  记录真值日志
%
%   M0Logs = append_truth_logs_m0(M0Logs, FrameState_t, Config)
%
%   功能：
%     1. 对所有 ordinary_target 记录真值到 TruthLogTarget（供误差统计）
%     2. 对所有活跃源记录到 TruthLogAll（供全局分析）

    B = Config.m0.num_bands;

    for b = 1:B
        apb = FrameState_t.active_per_band(b);

        if ~apb.has_source
            continue;
        end

        % --- 记录到 TruthLogAll ---
        entry_all.frame_id     = FrameState_t.frame_id;
        entry_all.band_id      = b;
        entry_all.source_type  = apb.source_type;
        entry_all.instance_id  = apb.instance_id;
        entry_all.template_key = apb.template_key;
        entry_all.true_pos_xy  = apb.true_pos_xy;
        entry_all.tx_power_dBm = apb.tx_power_dBm;

        if isempty(M0Logs.TruthLogAll)
            M0Logs.TruthLogAll = entry_all;
        else
            M0Logs.TruthLogAll(end+1) = entry_all;
        end

        % --- 对 ordinary_target 额外记录到 TruthLogTarget ---
        if apb.is_localization_target
            entry_tgt.frame_id     = FrameState_t.frame_id;
            entry_tgt.band_id      = b;
            entry_tgt.instance_id  = apb.instance_id;
            entry_tgt.template_key = apb.template_key;
            entry_tgt.true_pos_xy  = apb.true_pos_xy;
            entry_tgt.tx_power_dBm = apb.tx_power_dBm;

            if isempty(M0Logs.TruthLogTarget)
                M0Logs.TruthLogTarget = entry_tgt;
            else
                M0Logs.TruthLogTarget(end+1) = entry_tgt;
            end
        end
    end
end
