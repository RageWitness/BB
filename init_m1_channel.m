function [Bands, ChannelState, Config] = init_m1_channel(Config, APs)
% INIT_M1_CHANNEL  M1 信道模块总初始化入口
%
%   [Bands, ChannelState, Config] = init_m1_channel(Config, APs)
%
%   功能：
%     1. 填充 M1 默认配置（若缺失）
%     2. 生成 Bands 结构体
%     3. 初始化 ChannelState
%
%   输出：
%       Bands        - 频带参数结构体
%       ChannelState - 信道状态（初始为零）
%       Config       - 补充后的完整配置

    %% 1. 填充 M1 默认配置
    Config = fill_m1_defaults(Config);

    %% 2. 生成 Bands
    B = Config.m0.num_bands;
    Bands.B      = B;
    Bands.fc_Hz  = Config.m1.bands.fc_Hz;
    Bands.bw_Hz  = Config.m1.bands.bw_Hz;
    Bands.name   = Config.m1. bands.name;
    Bands.model  = Config.m1.channel.model_per_band;

    fprintf('[M1] 频带配置:\n');
    for b = 1:B
        fprintf('  Band %d: %s, fc=%.1f MHz, BW=%.1f kHz, model=%s\n', ...
            b, Bands.name{b}, Bands.fc_Hz(b)/1e6, Bands.bw_Hz(b)/1e3, Bands.model{b});
    end

    %% 3. 初始化 ChannelState
    ChannelState = init_m1_channel_state(Config, APs, Bands);

    fprintf('[M1] ====== M1 信道初始化完成 ======\n');
end


function Config = fill_m1_defaults(Config)
% FILL_M1_DEFAULTS  填充 M1 缺失的默认配置

    B = Config.m0.num_bands;

    % --- 频带定义 ---
    if ~isfield(Config, 'm1') || ~isfield(Config.m1, 'bands')
        Config.m1.bands.fc_Hz = [96e6, 89e6, 2.4e9, 2593e6];
        Config.m1.bands.bw_Hz = [200e3, 200e3, 40e6, 194e6];
        Config.m1.bands.name  = {'VHF-96M', 'VHF-89M', 'ISM-2.4G', 'S-2.6G'};
    end

    % --- 信道模型选择（统一 MWM 模型）---
    if ~isfield(Config.m1, 'channel') || ~isfield(Config.m1.channel, 'model_per_band')
        Config.m1.channel.model_per_band = {'mwm', 'mwm', 'mwm', 'mwm'};
    end

    % --- 建筑布局（新统一模型使用）---
    if ~isfield(Config.m1.channel, 'buildings') || isempty(Config.m1.channel.buildings)
        Config.m1.channel.buildings = get_default_buildings();
    end

    % --- legacy 保留字段（不再使用，仅防止下游读取报错）---
    % 旧 lognormal / COST231-WI 配置块已删除 (统一 MWM 模型不再需要)

    % --- 慢变 ---
    if ~isfield(Config.m1.channel, 'enableSlowDrift')
        Config.m1.channel.enableSlowDrift = false;
    end
    if ~isfield(Config.m1, 'slow')
        Config.m1.slow.global.rho   = 0.99;
        Config.m1.slow.global.sigma = 0.3;
        Config.m1.slow.ap.rho       = 0.995;
        Config.m1.slow.ap.sigma     = 0.2;
        Config.m1.slow.shadow.rho   = 0.98;
        Config.m1.slow.shadow.sigma = 4 * ones(1, B);
    end

    % --- 噪声 ---
    if ~isfield(Config.m1, 'noise')
        Config.m1.noise.mode           = 'psd_dBmHz';
        Config.m1.noise.n0_mode        = 'fixed';
        Config.m1.noise.n0_dBmHz       = -174 * ones(1, B);
        Config.m1.noise.noise_power_dBm = [];  % power_dBm 模式备用
        Config.m1.noise.sliding_rho    = 0.999;
        Config.m1.noise.sliding_sigma  = 0.5;
    end
end
