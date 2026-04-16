function [SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(Config_override)
% INIT_M0_NONOVERLAP  M0 模块总初始化入口
%
%   [SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap()
%   [SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(Config_override)
%
%   功能：
%     1. 生成默认配置（或合并用户覆盖）
%     2. 生成 AP 布局
%     3. 构建有效网格
%     4. 构建三类信源模板库
%     5. 初始化运行时状态
%     6. 初始化真值日志

    %% ===== 1. 默认配置 =====
    Config = default_config_m0();

    % 如果提供了覆盖参数，递归合并
    if nargin >= 1 && ~isempty(Config_override)
        Config = merge_struct(Config, Config_override);
    end

    %% ===== 2. AP 布局 =====
    APs = generate_aps(Config);

    %% ===== 3. 有效网格 =====
    GridValid = build_valid_grid_nonoverlap(Config, APs);

    %% ===== 4. 模板库 =====
    SourceTemplates = build_source_templates_nonoverlap(Config, GridValid);

    %% ===== 5. 运行时状态 =====
    M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates);

    %% ===== 6. 真值日志 =====
    M0Logs.TruthLogTarget = struct( ...
        'frame_id',     {}, ...
        'band_id',      {}, ...
        'instance_id',  {}, ...
        'template_key', {}, ...
        'true_pos_xy',  {}, ...
        'tx_power_dBm', {});

    M0Logs.TruthLogAll = struct( ...
        'frame_id',     {}, ...
        'band_id',      {}, ...
        'source_type',  {}, ...
        'instance_id',  {}, ...
        'template_key', {}, ...
        'true_pos_xy',  {}, ...
        'tx_power_dBm', {});

    fprintf('[M0] ====== M0 初始化完成 ======\n');
end


%% ===================== 内部函数 =====================

function Config = default_config_m0()
% DEFAULT_CONFIG_M0  M0 默认配置

    % --- 场景 ---
    Config.area.x_range = [0, 300];   % m
    Config.area.y_range = [0, 300];   % m

    % --- AP ---
    Config.ap.num_x = 4;
    Config.ap.num_y = 4;

    % --- 指纹网格 ---
    Config.fp.grid_step           = 3;    % m
    Config.fp.ap_exclusion_radius = 1.5;  % m

    % --- M0 基本参数 ---
    Config.m0.num_bands = 4;
    Config.m0.dt        = 1;        % 帧间隔 (s)
    Config.m0.T_total   = 500;      % 总帧数

    Config.m0.priority_order = {'trusted_fixed', 'prior_pos_known', ...
                                'prior_time_known', 'ordinary_target'};

    % --- trusted_fixed ---
    Config.m0.trusted.TrustedNum     = 2;
    Config.m0.trusted.lambda_arrival = 0.005;   % 到达率 (1/s)
    Config.m0.trusted.life_mode      = 'geom';  % 'geom' 或 'exp'
    Config.m0.trusted.life_param     = 0.05;    % geom: p_end; exp: mu(s)

    % --- prior ---
    Config.m0.prior.Sprior           = [2, 2, 1, 1];  % 每频带 prior_pos_known 数量
    Config.m0.prior.Tprior           = [1, 1, 1, 1];  % 每频带 prior_time_known 数量
    Config.m0.prior.schedule_mode    = 'periodic';     % 'periodic' 或 'table'
    Config.m0.prior.period_frames    = [50, 60, 70, 80];
    Config.m0.prior.duration_frames  = [10, 10, 8, 8];
    Config.m0.prior.phase_frames     = [0, 5, 10, 15];

    % --- ordinary_target ---
    Config.m0.target.Nocoop          = [3, 3, 2, 2];   % 每频带模板数
    Config.m0.target.lambda_arrival  = [0.02, 0.02, 0.015, 0.015];
    Config.m0.target.life_mode       = 'geom';          % 'geom' 或 'uniform'
    Config.m0.target.life_param      = 0.08;            % geom: p_end
    Config.m0.target.life_range      = [5, 30];         % uniform 模式下 [Lmin, Lmax]
    Config.m0.target.power_range_dBm = [50, 65;   ... % band 1
                                        50, 65;   ... % band 2
                                        50, 65;   ... % band 3
                                        50, 65];      % band 4
    Config.m0.target.position_mode   = 'uniform';       % 'uniform' 或 'hotspot'

    % --- M3 grouping 配置 ---
    Config.m3.grouping.enable = true;
    Config.m3.grouping.max_gap_frames = 2;

    % --- M3 相对功率特征配置（辅助证据） ---
    Config.m3.relative_power.enable = true;
    Config.m3.relative_power.topk = 3;
    Config.m3.relative_power.tau_valid_ap_dB = 3;

    % --- M3 trusted hard gate ---
    Config.m3.trusted.min_active_bands = 3;
    Config.m3.trusted.min_power_z = 1.5;
    Config.m3.trusted.max_power_jitter_dB = 3;
    Config.m3.trusted.min_position_score = 0.70;
    Config.m3.trusted.max_position_spread_m = 10;

    % --- M3 prior gate ---
    Config.m3.prior_pos.min_position_score = 0.60;
    Config.m3.prior_pos.min_margin = 0.10;
    Config.m3.prior_time.min_time_score = 0.60;
    Config.m3.prior_time.min_margin = 0.10;

    % --- M3 ordinary fallback ---
    Config.m3.ordinary.min_detect_z = 0.8;
    Config.m3.ordinary.min_duration_frames = 3;
    Config.m3.ordinary.min_valid_ap = 2;

    % --- M3 route ---
    Config.m3.route.hold_enable = true;

    % --- M4 掩码 shape 配置 ---
    Config.m4.distance_mode = 'shape_scale_masked'; % 可选: shape_scale/shape_scale_prob/shape_scale_hn
    Config.m4.masked_shape.enable = true;
    Config.m4.masked_shape.tau_low_dB = 3;
    Config.m4.masked_shape.tau_high_dB = 10;
    Config.m4.masked_shape.eps = 1e-12;
    Config.m4.lambda_shape = 0.7;
    Config.m4.lambda_resid = 0.3;
end


function APs = generate_aps(Config)
% GENERATE_APS  生成均匀分布的 AP
    nx = Config.ap.num_x;
    ny = Config.ap.num_y;
    x_range = Config.area.x_range;
    y_range = Config.area.y_range;

    margin_x = (x_range(2) - x_range(1)) / (2*nx);
    margin_y = (y_range(2) - y_range(1)) / (2*ny);

    ax = linspace(x_range(1) + margin_x, x_range(2) - margin_x, nx);
    ay = linspace(y_range(1) + margin_y, y_range(2) - margin_y, ny);

    [AX, AY] = meshgrid(ax, ay);
    APs.pos_xy = [AX(:), AY(:)];
    APs.num    = nx * ny;

    fprintf('[M0] AP 数量: %d (%dx%d 均匀分布)\n', APs.num, nx, ny);
end


function S = merge_struct(S, T)
% MERGE_STRUCT  递归合并结构体 T 的字段覆盖到 S
    fnames = fieldnames(T);
    for i = 1:numel(fnames)
        fn = fnames{i};
        if isfield(S, fn) && isstruct(S.(fn)) && isstruct(T.(fn))
            S.(fn) = merge_struct(S.(fn), T.(fn));
        else
            S.(fn) = T.(fn);
        end
    end
end
