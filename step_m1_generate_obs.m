function [ChannelState, ObsFrame] = step_m1_generate_obs( ...
    FrameState_t, APs, Bands, ChannelState, Config, t)
% STEP_M1_GENERATE_OBS  M1 单帧总入口：生成当前帧 RSS 观测
%
%   [ChannelState, ObsFrame] = step_m1_generate_obs(...)
%
%   流程：
%     1. update_m1_slow_drift    — 更新慢变状态（若开启）
%     2. update_m1_noise_state   — 更新噪声 n0（若滑动）
%     3. 对 b=1:B 调用 generate_band_obs_m1
%     4. assemble_obs_frame_m1   — 组装输出
%
%   输入：
%       FrameState_t  - M0 当前帧真值状态
%       APs           - AP 信息
%       Bands         - 频带参数
%       ChannelState  - 信道状态
%       Config        - 配置
%       t             - 当前帧号
%
%   输出：
%       ChannelState  - 更新后的信道状态
%       ObsFrame      - 当前帧观测

    B = Bands.B;

    %% 1. 更新慢变状态
    ChannelState = update_m1_slow_drift(ChannelState, Config);

    %% 2. 更新噪声状态
    ChannelState = update_m1_noise_state(ChannelState, Config, t);

    %% 3. 生成各频带观测
    BandObsList = cell(1, B);
    for b = 1:B
        BandObsList{b} = generate_band_obs_m1( ...
            FrameState_t.active_per_band(b), APs, Bands, ...
            ChannelState, Config, b);
    end

    %% 4. 组装 ObsFrame
    ObsFrame = assemble_obs_frame_m1(BandObsList, FrameState_t, Config);
end
