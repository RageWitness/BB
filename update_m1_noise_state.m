function ChannelState = update_m1_noise_state(ChannelState, Config, t)
% UPDATE_M1_NOISE_STATE  更新噪声 n0 状态
%
%   ChannelState = update_m1_noise_state(ChannelState, Config, t)
%
%   若 n0_mode = 'fixed'：不做操作
%   若 n0_mode = 'sliding'：AR(1) 更新 noise_n0_dBmHz

    if strcmp(Config.m1.noise.n0_mode, 'fixed')
        return;
    end

    % sliding 模式：AR(1)
    rho   = Config.m1.noise.sliding_rho;
    sigma = Config.m1.noise.sliding_sigma;
    B = numel(ChannelState.noise_n0_dBmHz);

    ChannelState.noise_n0_dBmHz = rho * ChannelState.noise_n0_dBmHz + ...
                                  (1 - rho) * Config.m1.noise.n0_dBmHz + ...
                                  sigma * randn(1, B);
end
