function M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates)
% INIT_RUNTIME_STATE_M0_NONOVERLAP  初始化 M0 运行时状态（四类源）

    B = Config.m0.num_bands;

    M0State.frame_id = 0;
    M0State.time_sec = 0;
    M0State.next_instance_id = 1;

    % 实例池（轻量化记账，仅供生命周期追踪）
    M0State.instances = struct( ...
        'instance_id',    {}, ...
        'template_key',   {}, ...
        'source_type',    {}, ...
        'band_list',      {}, ...
        'true_position',  {}, ...
        'true_tx_power',  {}, ...
        'tx_power_by_band', {}, ...
        'location_prior', {}, ...
        'power_prior',    {}, ...
        'life_remaining', {}, ...
        'is_active',      {});

    % 每频带活跃源
    empty_band.has_source   = false;
    empty_band.instance_id  = 0;
    empty_band.template_key = '';
    for b = 1:B
        M0State.active_per_band(b) = empty_band;
    end

    % broadband_cal 跟踪
    BC = Config.m0.source.broadband_cal;
    if isfield(BC, 'schedule_mode') && strcmp(BC.schedule_mode, 'manual')
        BC_count = numel(BC.manual_schedule);
    else
        BC_count = BC.count;
    end
    M0State.broadband_active  = false(1, BC_count);
    M0State.broadband_life    = zeros(1, BC_count);
    M0State.broadband_inst_id = zeros(1, BC_count);

    % target 每频带跟踪
    M0State.target_active_band = false(1, B);
    M0State.target_life_band   = zeros(1, B);
    M0State.target_inst_band   = zeros(1, B);

    % opportunistic 每频带跟踪（一频带至多一个机遇源实例）
    M0State.opp_active_band = false(1, B);
    M0State.opp_life_band   = zeros(1, B);
    M0State.opp_inst_band   = zeros(1, B);

    fprintf('[M0] 运行时状态初始化完成\n');
end
