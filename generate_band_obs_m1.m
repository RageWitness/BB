function BandObs = generate_band_obs_m1(active_band_state, APs, Bands, ...
    ChannelState, Config, band_id)
% GENERATE_BAND_OBS_M1  生成单频带当前帧的 RSS 观测
%
%   BandObs = generate_band_obs_m1(active_band_state, APs, Bands,
%                                   ChannelState, Config, band_id)
%
%   若无源：输出噪声底
%   若有源：信号功率 + 接收机噪声底功率
%
%   输出：
%       BandObs.Y_dBm              - (M x 1) RSS 观测 (dBm)
%       BandObs.Y_lin              - (M x 1) RSS 观测 (线性 mW)
%       BandObs.meta.has_source
%       BandObs.meta.pathloss_dB   - (M x 1) 路径损耗
%       BandObs.meta.rx_signal_dBm - (M x 1) 纯信号接收功率
%       BandObs.meta.noise_power_dBm

    b = band_id;
    M = APs.num;

    % --- 计算噪声功率 ---
    if strcmp(Config.m1.noise.mode, 'psd_dBmHz')
        N_dBm = ChannelState.noise_n0_dBmHz(b) + 10 * log10(Bands.bw_Hz(b));
    else
        N_dBm = Config.m1.noise.noise_power_dBm(b);
    end
    N_lin = 10^(N_dBm / 10);  % mW

    BandObs.meta.noise_power_dBm = N_dBm;

    if ~active_band_state.has_source
        % ===== 无源：纯噪声底 =====
        BandObs.Y_lin  = N_lin * ones(M, 1);
        BandObs.Y_dBm  = 10 * log10(BandObs.Y_lin);
        BandObs.meta.has_source    = false;
        BandObs.meta.pathloss_dB   = zeros(M, 1);
        BandObs.meta.rx_signal_dBm = -inf * ones(M, 1);
        BandObs.meta.noise_model   = 'mean_power_floor';
        return;
    end

    % ===== 有源 =====
    src_xy   = active_band_state.true_pos_xy;
    tx_power = active_band_state.tx_power_dBm;

    % 1. 计算源到各 AP 距离（仅供 meta 记录）
    dist_m = sqrt(sum((APs.pos_xy - src_xy).^2, 2));  %#ok<NASGU>

    % 2. 根据频带模型计算路径损耗
    band_model = Config.m1.channel.model_per_band{b};
    if strcmp(band_model, 'lognormal')
        [PL_dB, ~] = compute_pathloss_lognormal(src_xy, APs.pos_xy, Bands.fc_Hz(b), Config);
    else
        buildings = Config.m1.channel.buildings;
        [PL_dB, ~] = compute_pathloss_mwm(src_xy, APs.pos_xy, Bands.fc_Hz(b), buildings);
    end

    % 3. 叠加慢变偏置（MWM 模型下阴影由 sigma 残差承担，此处仅全局/AP 偏置）
    drift_dB = zeros(M, 1);
    if Config.m1.channel.enableSlowDrift
        drift_dB = drift_dB + ChannelState.global_bias_dB(b);
        drift_dB = drift_dB + ChannelState.ap_bias_dB(:, b);
        drift_dB = drift_dB + ChannelState.shadow_state_dB(:, b);
    end

    % 4. 接收信号功率
    rx_signal_dBm = tx_power - PL_dB + drift_dB;  % M x 1
    S_lin = 10.^(rx_signal_dBm / 10);              % mW

    % 5. Add receiver noise power floor in linear power domain.
    % N_lin is already a power in mW.
    Y_lin = S_lin + N_lin;

    % 确保非负（物理功率不能为负）
    Y_lin = max(Y_lin, 1e-20);

    % 6. 转回 dBm
    Y_dBm = 10 * log10(Y_lin);

    % 组装输出
    BandObs.Y_dBm  = Y_dBm;
    BandObs.Y_lin  = Y_lin;
    BandObs.meta.has_source    = true;
    BandObs.meta.pathloss_dB   = PL_dB;
    BandObs.meta.rx_signal_dBm = rx_signal_dBm;
    BandObs.meta.noise_model   = 'mean_power_floor';
end
