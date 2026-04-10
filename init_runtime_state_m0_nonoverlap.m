function M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates)
% INIT_RUNTIME_STATE_M0_NONOVERLAP  初始化 M0 运行时状态
%
%   M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates)
%
%   输出：
%       M0State.frame_id          - 当前帧号 (初始 0)
%       M0State.time_sec          - 当前时间 (s)
%       M0State.next_instance_id  - 下一个实例 ID 计数器
%       M0State.instances         - 实例数组 (初始为空)
%       M0State.active_per_band   - (1 x B) 每频带活跃源信息
%       M0State.trusted_active    - (1 x TrustedNum) 各 trusted 模板的活跃状态

    B = Config.m0.num_bands;

    M0State.frame_id = 0;
    M0State.time_sec = 0;
    M0State.next_instance_id = 1;

    % 实例池 —— 使用结构体数组，初始为空
    M0State.instances = struct( ...
        'instance_id',          {}, ...
        'template_key',         {}, ...
        'source_type',          {}, ...
        'source_subtype',       {}, ...
        'band_id_list',         {}, ...
        'true_pos_xy',          {}, ...
        'tx_power_dBm',         {}, ...
        'life_remaining',       {}, ...
        'is_active',            {}, ...
        'system_knows_position',{}, ...
        'system_knows_time',    {});

    % 每频带活跃源 —— 初始均为空
    empty_band.has_source   = false;
    empty_band.instance_id  = 0;
    empty_band.template_key = '';
    for b = 1:B
        M0State.active_per_band(b) = empty_band;
    end

    % trusted 模板活跃状态追踪
    TrustedNum = Config.m0.trusted.TrustedNum;
    M0State.trusted_active = false(1, TrustedNum);
    M0State.trusted_life   = zeros(1, TrustedNum);

    % ordinary_target 每频带是否有持续中的目标
    M0State.target_active_band = false(1, B);
    M0State.target_life_band   = zeros(1, B);
    M0State.target_inst_band   = zeros(1, B);  % 实例 ID

    fprintf('[M0] 运行时状态初始化完成\n');
end
