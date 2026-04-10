function ObsFrame = assemble_obs_frame_m1(BandObsList, FrameState_t, Config)
% ASSEMBLE_OBS_FRAME_M1  将各频带观测结果组装成 ObsFrame
%
%   ObsFrame = assemble_obs_frame_m1(BandObsList, FrameState_t, Config)
%
%   输入：
%       BandObsList - (1 x B) cell 数组，每个元素为 BandObs 结构体
%       FrameState_t - M0 输出的帧状态
%       Config       - 配置
%
%   输出：
%       ObsFrame.frame_id
%       ObsFrame.time_sec
%       ObsFrame.Y_dBm     - (M x B) RSS 观测矩阵 (dBm)
%       ObsFrame.Y_lin     - (M x B) RSS 观测矩阵 (线性 mW)
%       ObsFrame.meta      - 各频带元信息

    B = Config.m0.num_bands;
    M = numel(BandObsList{1}.Y_dBm);

    ObsFrame.frame_id = FrameState_t.frame_id;
    ObsFrame.time_sec = FrameState_t.time_sec;

    ObsFrame.Y_dBm = zeros(M, B);
    ObsFrame.Y_lin = zeros(M, B);

    for b = 1:B
        bo = BandObsList{b};
        ObsFrame.Y_dBm(:, b) = bo.Y_dBm;
        ObsFrame.Y_lin(:, b) = bo.Y_lin;

        % meta：来自 BandObs + FrameState
        apb = FrameState_t.active_per_band(b);
        ObsFrame.meta.active_per_band(b).has_source       = apb.has_source;
        ObsFrame.meta.active_per_band(b).template_key     = apb.template_key;
        ObsFrame.meta.active_per_band(b).source_type      = apb.source_type;
        ObsFrame.meta.active_per_band(b).true_pos_xy      = apb.true_pos_xy;
        ObsFrame.meta.active_per_band(b).tx_power_dBm     = apb.tx_power_dBm;
        ObsFrame.meta.active_per_band(b).pathloss_dB_per_ap   = bo.meta.pathloss_dB;
        ObsFrame.meta.active_per_band(b).rx_signal_dBm_per_ap = bo.meta.rx_signal_dBm;
        ObsFrame.meta.active_per_band(b).noise_power_dBm      = bo.meta.noise_power_dBm;
    end
end
