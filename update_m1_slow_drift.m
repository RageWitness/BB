function ChannelState = update_m1_slow_drift(ChannelState, Config)
% UPDATE_M1_SLOW_DRIFT  更新慢变状态（AR(1) 过程）
%
%   ChannelState = update_m1_slow_drift(ChannelState, Config)
%
%   若 enableSlowDrift = false，不做任何操作。
%   若 enableSlowDrift = true，按 AR(1) 更新三层慢变：
%     1. global_bias_dB(b)  — 频带全局漂移
%     2. ap_bias_dB(m,b)    — AP 端偏置
%     3. shadow_state_dB(m,b) — 阴影慢变

    if ~Config.m1.channel.enableSlowDrift
        return;
    end

    [M, B] = size(ChannelState.ap_bias_dB);

    % 1. 频带全局慢变
    rho_g   = Config.m1.slow.global.rho;
    sigma_g = Config.m1.slow.global.sigma;
    ChannelState.global_bias_dB = rho_g * ChannelState.global_bias_dB + ...
                                  sigma_g * randn(1, B);

    % 2. AP 端慢变
    rho_a   = Config.m1.slow.ap.rho;
    sigma_a = Config.m1.slow.ap.sigma;
    ChannelState.ap_bias_dB = rho_a * ChannelState.ap_bias_dB + ...
                              sigma_a * randn(M, B);

    % 3. 阴影慢变
    rho_s   = Config.m1.slow.shadow.rho;
    sigma_s = Config.m1.slow.shadow.sigma;  % 1 x B
    for b = 1:B
        ChannelState.shadow_state_dB(:, b) = ...
            rho_s * ChannelState.shadow_state_dB(:, b) + ...
            sigma_s(b) * sqrt(1 - rho_s^2) * randn(M, 1);
    end
end
