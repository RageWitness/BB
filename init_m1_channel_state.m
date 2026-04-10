 function ChannelState = init_m1_channel_state(Config, APs, Bands)
% INIT_M1_CHANNEL_STATE  初始化 M1 信道状态
%
%   ChannelState = init_m1_channel_state(Config, APs, Bands)
%
%   输出：
%       ChannelState.global_bias_dB   - (1 x B) 频带全局慢变偏置
%       ChannelState.ap_bias_dB       - (M x B) AP 端慢变偏置
%       ChannelState.shadow_state_dB  - (M x B) 阴影慢变状态
%       ChannelState.noise_n0_dBmHz   - (1 x B) 当前噪声 PSD

    M = APs.num;
    B = Bands.B;

    ChannelState.global_bias_dB  = zeros(1, B);
    ChannelState.ap_bias_dB      = zeros(M, B);
    ChannelState.shadow_state_dB = zeros(M, B);

    % 噪声 PSD 初始化
    if strcmp(Config.m1.noise.mode, 'psd_dBmHz')
        ChannelState.noise_n0_dBmHz = Config.m1.noise.n0_dBmHz;
    else
        % power_dBm 模式下不需要 PSD，但仍初始化为备用
        ChannelState.noise_n0_dBmHz = -174 * ones(1, B);
    end

    fprintf('[M1] ChannelState 初始化完成 (M=%d, B=%d)\n', M, B);
end
