function SourceTemplates = build_source_templates_nonoverlap(Config, GridValid)
% BUILD_SOURCE_TEMPLATES_NONOVERLAP  构建三类信源模板库
%
%   SourceTemplates = build_source_templates_nonoverlap(Config, GridValid)
%
%   模板类别：
%     trusted(j)         — 宽带可信标校源
%     prior_pos(b,s)     — 位置已知先验源
%     prior_time(b,q)    — 时间已知先验源
%     target(b,n)        — 普通待定位目标
%
%   输入：
%       Config    - 配置结构体
%       GridValid - 有效网格结构体
%
%   输出：
%       SourceTemplates - 模板库结构体

    B = Config.m0.num_bands;

    %% ===== A. trusted_fixed =====
    TrustedNum = Config.m0.trusted.TrustedNum;
    for j = 1:TrustedNum
        tpl.source_type     = 'trusted_fixed';
        tpl.template_key    = sprintf('anchor_%03d', j);
        tpl.bands_covered   = 1:B;
        % 在场景中均匀选取固定位置
        idx = randi(GridValid.Nvalid);
        tpl.fixed_pos_xy    = GridValid.xy(idx, :);
        tpl.tx_power_dBm    = 15;
        tpl.lambda_arrival  = Config.m0.trusted.lambda_arrival;
        tpl.life_mode       = Config.m0.trusted.life_mode;
        tpl.life_param      = Config.m0.trusted.life_param;
        tpl.credibility_prior = 0.99;
        SourceTemplates.trusted(j) = tpl;
    end

    %% ===== B1. prior_pos_known =====
    for b = 1:B
        Sprior_b = Config.m0.prior.Sprior(b);
        for s = 1:Sprior_b
            tpl.source_type     = 'prior';
            tpl.source_subtype  = 'prior_pos_known';
            tpl.template_key    = sprintf('prior_pos_b%d_%03d', b, s);
            tpl.band_id         = b;
            idx = randi(GridValid.Nvalid);
            tpl.fixed_pos_xy    = GridValid.xy(idx, :);
            tpl.tx_power_dBm    = 0;
            tpl.schedule_mode   = Config.m0.prior.schedule_mode;

            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.period_frames   = Config.m0.prior.period_frames(b);
                tpl.duration_frames = Config.m0.prior.duration_frames(b);
                tpl.phase_frames    = Config.m0.prior.phase_frames(b) + (s-1)*3;
            else
                % timetable 模式：生成随机出现时间集合
                T_total = Config.m0.T_total;
                n_events = max(1, round(T_total / Config.m0.prior.period_frames(b)));
                schedule_set = sort(randperm(T_total, min(n_events * Config.m0.prior.duration_frames(b), T_total)));
                tpl.schedule_set = schedule_set;
            end

            SourceTemplates.prior_pos(b,s) = tpl;
        end
    end

    %% ===== B2. prior_time_known =====
    for b = 1:B
        Tprior_b = Config.m0.prior.Tprior(b);
        for q = 1:Tprior_b
            tpl.source_type     = 'prior';
            tpl.source_subtype  = 'prior_time_known';
            tpl.template_key    = sprintf('prior_time_b%d_%03d', b, q);
            tpl.band_id         = b;
            tpl.tx_power_dBm    = 0;
            tpl.location_mode   = 'unknown_to_system';
            tpl.schedule_mode   = Config.m0.prior.schedule_mode;

            if strcmp(tpl.schedule_mode, 'periodic')
                tpl.period_frames   = Config.m0.prior.period_frames(b);
                tpl.duration_frames = Config.m0.prior.duration_frames(b);
                % 与 prior_pos 错开相位
                tpl.phase_frames    = Config.m0.prior.phase_frames(b) + ...
                                      Config.m0.prior.Sprior(b)*3 + (q-1)*3;
            else
                T_total = Config.m0.T_total;
                n_events = max(1, round(T_total / Config.m0.prior.period_frames(b)));
                schedule_set = sort(randperm(T_total, min(n_events * Config.m0.prior.duration_frames(b), T_total)));
                tpl.schedule_set = schedule_set;
            end

            SourceTemplates.prior_time(b,q) = tpl;
        end
    end

    %% ===== C. ordinary_target =====
    for b = 1:B
        Nocoop_b = Config.m0.target.Nocoop(b);
        for n = 1:Nocoop_b
            tpl.source_type      = 'ordinary_target';
            tpl.source_subtype   = 'ordinary_target';
            tpl.template_key     = sprintf('target_b%d_%03d', b, n);
            tpl.band_id          = b;
            tpl.lambda_arrival   = Config.m0.target.lambda_arrival(b);
            tpl.life_mode        = Config.m0.target.life_mode;
            tpl.life_param       = Config.m0.target.life_param;
            tpl.tx_power_range_dBm = Config.m0.target.power_range_dBm(b, :);
            tpl.position_mode    = Config.m0.target.position_mode;
            SourceTemplates.target(b,n) = tpl;
        end
    end

    fprintf('[M0] 模板库构建完成: trusted=%d, prior_pos=%s, prior_time=%s, target=%s\n', ...
        TrustedNum, ...
        mat2str(Config.m0.prior.Sprior), ...
        mat2str(Config.m0.prior.Tprior), ...
        mat2str(Config.m0.target.Nocoop));
end
