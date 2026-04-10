function M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates)
% INIT_RUNTIME_STATE_M0_NONOVERLAP  初始化 M0 运行时状态（含先验对象字段）
%
%   M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates)

    B = Config.m0.num_bands;

    M0State.frame_id = 0;
    M0State.time_sec = 0;
    M0State.next_instance_id = 1;

    % 空的先验对象模板（用于初始化结构体字段）
    empty_pos_prior.mode                = '';
    empty_pos_prior.candidate_positions = [];
    empty_pos_prior.match_radius        = [];

    empty_time_prior.mode            = '';
    empty_time_prior.schedule        = [];
    empty_time_prior.period_frames   = [];
    empty_time_prior.duration_frames = [];
    empty_time_prior.phase_frames    = [];

    % 实例池 —— 结构体数组，包含先验对象和模板特征
    M0State.instances = struct( ...
        'instance_id',           {}, ...
        'template_key',          {}, ...
        'source_type',           {}, ...
        'source_subtype',        {}, ...
        'band_id_list',          {}, ...
        'true_pos_xy',           {}, ...
        'tx_power_dBm',          {}, ...
        'life_remaining',        {}, ...
        'is_active',             {}, ...
        'position_prior',        {}, ...
        'time_prior',            {}, ...
        'power_nominal_dBm',     {}, ...
        'power_range_dBm',       {}, ...
        'power_stability_level', {}, ...
        'credibility_prior',     {}, ...
        'upgrade_potential',     {});

    % 每频带活跃源 —— 初始均为空
    empty_band.has_source   = false;
    empty_band.instance_id  = 0;
    empty_band.template_key = '';
    for b = 1:B
        M0State.active_per_band(b) = empty_band;
    end

    % trusted 模板活跃状态追踪
    TrustedNum = Config.m0.trusted.TrustedNum;
    M0State.trusted_active  = false(1, TrustedNum);
    M0State.trusted_life    = zeros(1, TrustedNum);
    M0State.trusted_inst_id = zeros(1, TrustedNum);

    % ordinary_target 每频带追踪
    M0State.target_active_band = false(1, B);
    M0State.target_life_band   = zeros(1, B);
    M0State.target_inst_band   = zeros(1, B);

    fprintf('[M0] 运行时状态初始化完成\n');
end
