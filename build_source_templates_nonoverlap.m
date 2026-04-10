function SourceTemplates = build_source_templates_nonoverlap(Config, GridValid)
% BUILD_SOURCE_TEMPLATES_NONOVERLAP  构建三类信源模板库（含先验对象）
%
%   SourceTemplates = build_source_templates_nonoverlap(Config, GridValid)
%
%   每个模板包含两层信息：
%     A. 静态模板特征：source_type, band_mask, power, credibility 等
%     B. 先验对象：position_prior, time_prior
%
%   模板类别：
%     trusted(j)         — 宽带可信标校源
%     prior_pos(b,s)     — 位置已知先验源
%     prior_time(b,q)    — 时间已知先验源
%     target(b,n)        — 普通待定位目标

    B = Config.m0.num_bands;

    %% ===== A. trusted_fixed =====
    TrustedNum = Config.m0.trusted.TrustedNum;
    for j = 1:TrustedNum
        tpl = struct();

        % --- 静态模板特征 ---
        tpl.source_type          = 'trusted_fixed';
        tpl.source_subtype       = 'trusted_fixed';
        tpl.template_namespace   = 'trusted_fixed';
        tpl.template_key         = sprintf('anchor_%03d', j);
        tpl.band_mask            = ones(1, B);         % 覆盖所有频带
        tpl.bands_covered        = 1:B;

        % 功率
        tpl.power_nominal_dBm    = 15;
        tpl.power_range_dBm      = [15, 15];           % 固定功率
        tpl.power_stability_level = 3;                  % 最稳定
        tpl.tx_power_dBm         = 15;

        % 到达与寿命
        tpl.lambda_arrival       = Config.m0.trusted.lambda_arrival;
        tpl.life_mode            = Config.m0.trusted.life_mode;
        tpl.life_param           = Config.m0.trusted.life_param;
        tpl.time_pattern_type    = 'poisson';
        tpl.expected_duration_range = [10, 50];         % 预期持续帧数范围

        % 位置（固定）
        idx = randi(GridValid.Nvalid);
        tpl.fixed_pos_xy         = GridValid.xy(idx, :);

        % --- 位置先验对象 ---
        tpl.position_prior.mode                = 'fixed_known';
        tpl.position_prior.candidate_positions = tpl.fixed_pos_xy;  % 1x2 精确坐标
        tpl.position_prior.match_radius        = Config.fp.grid_step;  % 匹配半径

        % --- 时间先验对象 ---
        tpl.time_prior.mode            = 'unknown';     % 到达时间系统不知道
        tpl.time_prior.schedule        = [];
        tpl.time_prior.period_frames   = [];
        tpl.time_prior.duration_frames = [];
        tpl.time_prior.phase_frames    = [];

        % 可信度与晋升潜力
        tpl.credibility_prior    = 0.99;
        tpl.upgrade_potential    = 'none';              % 已经是最高级

        SourceTemplates.trusted(j) = tpl;
    end

    %% ===== B1. prior_pos_known =====
    for b = 1:B
        Sprior_b = Config.m0.prior.Sprior(b);
        for s = 1:Sprior_b
            tpl = struct();

            % --- 静态模板特征 ---
            tpl.source_type          = 'prior';
            tpl.source_subtype       = 'prior_pos_known';
            tpl.template_namespace   = 'prior';
            tpl.template_key         = sprintf('prior_pos_b%d_%03d', b, s);
            tpl.band_mask            = zeros(1, B); tpl.band_mask(b) = 1;
            tpl.band_id              = b;

            % 功率
            tpl.power_nominal_dBm    = 0;
            tpl.power_range_dBm      = [0, 0];
            tpl.power_stability_level = 3;
            tpl.tx_power_dBm         = 0;

            % 时间模式
            tpl.time_pattern_type    = 'scheduled';
            tpl.schedule_mode        = Config.m0.prior.schedule_mode;
            tpl.expected_duration_range = [Config.m0.prior.duration_frames(b), ...
                                           Config.m0.prior.duration_frames(b)];

            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.period_frames   = Config.m0.prior.period_frames(b);
                tpl.duration_frames = Config.m0.prior.duration_frames(b);
                tpl.phase_frames    = Config.m0.prior.phase_frames(b) + (s-1)*3;
            else
                T_total = Config.m0.T_total;
                n_events = max(1, round(T_total / Config.m0.prior.period_frames(b)));
                tpl.schedule_set = sort(randperm(T_total, ...
                    min(n_events * Config.m0.prior.duration_frames(b), T_total)));
            end

            % 位置（固定已知）
            idx = randi(GridValid.Nvalid);
            tpl.fixed_pos_xy         = GridValid.xy(idx, :);

            % --- 位置先验对象 ---
            tpl.position_prior.mode                = 'fixed_known';
            tpl.position_prior.candidate_positions = tpl.fixed_pos_xy;
            tpl.position_prior.match_radius        = Config.fp.grid_step;

            % --- 时间先验对象 ---
            tpl.time_prior.mode = tpl.schedule_mode;    % 'periodic' 或 'table'
            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.time_prior.schedule        = [];
                tpl.time_prior.period_frames   = tpl.period_frames;
                tpl.time_prior.duration_frames = tpl.duration_frames;
                tpl.time_prior.phase_frames    = tpl.phase_frames;
            else
                tpl.time_prior.schedule        = tpl.schedule_set;
                tpl.time_prior.period_frames   = [];
                tpl.time_prior.duration_frames = [];
                tpl.time_prior.phase_frames    = [];
            end

            % 可信度与晋升潜力
            tpl.credibility_prior    = 0.75 + 0.15 * rand();  % 0.75~0.9
            tpl.upgrade_potential    = 'none';

            SourceTemplates.prior_pos(b,s) = tpl;
        end
    end

    %% ===== B2. prior_time_known =====
    for b = 1:B
        Tprior_b = Config.m0.prior.Tprior(b);
        for q = 1:Tprior_b
            tpl = struct();

            % --- 静态模板特征 ---
            tpl.source_type          = 'prior';
            tpl.source_subtype       = 'prior_time_known';
            tpl.template_namespace   = 'prior';
            tpl.template_key         = sprintf('prior_time_b%d_%03d', b, q);
            tpl.band_mask            = zeros(1, B); tpl.band_mask(b) = 1;
            tpl.band_id              = b;

            % 功率
            tpl.power_nominal_dBm    = 0;
            tpl.power_range_dBm      = [0, 0];
            tpl.power_stability_level = 2;
            tpl.tx_power_dBm         = 0;

            % 时间模式
            tpl.time_pattern_type    = 'scheduled';
            tpl.location_mode        = 'unknown_to_system';
            tpl.schedule_mode        = Config.m0.prior.schedule_mode;
            tpl.expected_duration_range = [Config.m0.prior.duration_frames(b), ...
                                           Config.m0.prior.duration_frames(b)];

            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.period_frames   = Config.m0.prior.period_frames(b);
                tpl.duration_frames = Config.m0.prior.duration_frames(b);
                tpl.phase_frames    = Config.m0.prior.phase_frames(b) + ...
                                      Config.m0.prior.Sprior(b)*3 + (q-1)*3;
            else
                T_total = Config.m0.T_total;
                n_events = max(1, round(T_total / Config.m0.prior.period_frames(b)));
                tpl.schedule_set = sort(randperm(T_total, ...
                    min(n_events * Config.m0.prior.duration_frames(b), T_total)));
            end

            % --- 位置先验对象 —— 系统不知道位置 ---
            tpl.position_prior.mode                = 'time_known_position_unknown';
            tpl.position_prior.candidate_positions = [];  % 无候选
            tpl.position_prior.match_radius        = [];

            % --- 时间先验对象 ---
            tpl.time_prior.mode = tpl.schedule_mode;
            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.time_prior.schedule        = [];
                tpl.time_prior.period_frames   = tpl.period_frames;
                tpl.time_prior.duration_frames = tpl.duration_frames;
                tpl.time_prior.phase_frames    = tpl.phase_frames;
            else
                tpl.time_prior.schedule        = tpl.schedule_set;
                tpl.time_prior.period_frames   = [];
                tpl.time_prior.duration_frames = [];
                tpl.time_prior.phase_frames    = [];
            end

            % 可信度与晋升潜力
            tpl.credibility_prior    = 0.6 + 0.2 * rand();  % 0.6~0.8
            tpl.upgrade_potential    = 'none';

            SourceTemplates.prior_time(b,q) = tpl;
        end
    end

    %% ===== C. ordinary_target =====
    for b = 1:B
        Nocoop_b = Config.m0.target.Nocoop(b);
        for n = 1:Nocoop_b
            tpl = struct();

            % --- 静态模板特征 ---
            tpl.source_type          = 'ordinary_target';
            tpl.source_subtype       = 'ordinary_target';
            tpl.template_namespace   = 'ordinary_target';
            tpl.template_key         = sprintf('target_b%d_%03d', b, n);
            tpl.band_mask            = zeros(1, B); tpl.band_mask(b) = 1;
            tpl.band_id              = b;

            % 功率
            tpl.power_nominal_dBm    = mean(Config.m0.target.power_range_dBm(b, :));
            tpl.power_range_dBm      = Config.m0.target.power_range_dBm(b, :);
            tpl.power_stability_level = 1;              % 最不稳定
            tpl.tx_power_range_dBm   = Config.m0.target.power_range_dBm(b, :);

            % 到达与寿命
            tpl.lambda_arrival       = Config.m0.target.lambda_arrival(b);
            tpl.life_mode            = Config.m0.target.life_mode;
            tpl.life_param           = Config.m0.target.life_param;
            tpl.time_pattern_type    = 'poisson';
            tpl.expected_duration_range = Config.m0.target.life_range;
            tpl.position_mode        = Config.m0.target.position_mode;

            % --- 位置先验对象 —— 完全未知 ---
            tpl.position_prior.mode                = 'unknown';
            tpl.position_prior.candidate_positions = [];
            tpl.position_prior.match_radius        = [];

            % --- 时间先验对象 —— 完全未知 ---
            tpl.time_prior.mode            = 'unknown';
            tpl.time_prior.schedule        = [];
            tpl.time_prior.period_frames   = [];
            tpl.time_prior.duration_frames = [];
            tpl.time_prior.phase_frames    = [];

            % 可信度与晋升潜力
            tpl.credibility_prior    = 0.1 + 0.2 * rand();  % 0.1~0.3
            tpl.upgrade_potential    = assign_upgrade_potential(b, n, Config);

            SourceTemplates.target(b,n) = tpl;
        end
    end

    fprintf('[M0] 模板库构建完成: trusted=%d, prior_pos=%s, prior_time=%s, target=%s\n', ...
        TrustedNum, ...
        mat2str(Config.m0.prior.Sprior), ...
        mat2str(Config.m0.prior.Tprior), ...
        mat2str(Config.m0.target.Nocoop));
end


function up = assign_upgrade_potential(b, n, Config)
% ASSIGN_UPGRADE_POTENTIAL  根据模板在频带内的序号分配晋升潜力
%   第 1 个模板潜力最高（功率可能偏大、重现性可能更强）
    Nocoop_b = Config.m0.target.Nocoop(b);
    ratio = n / Nocoop_b;
    if ratio <= 0.33
        up = 'high';
    elseif ratio <= 0.66
        up = 'medium';
    else
        up = 'low';
    end
end
