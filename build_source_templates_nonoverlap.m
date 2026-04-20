function SourceTemplates = build_source_templates_nonoverlap(Config, GridValid)
% BUILD_SOURCE_TEMPLATES_NONOVERLAP  构建四类统一信源模板库
%
%   四类源（优先级从高到低）：
%     broadband_cal  — 宽带标校源（全频带）
%     opportunistic — 机遇源（模板在运行时动态生成实例与先验）
%     target        — 待定位源
%     persistent_cal — 持续存在的标校源（always-on 保底）

    B = Config.m0.num_bands;

    %% ===== 1. broadband_cal（宽带标校源，覆盖所有频带） =====
    BC = Config.m0.source.broadband_cal;
    for j = 1:BC.count
        tpl = struct();
        tpl.source_type        = 'broadband_cal';
        tpl.source_subtype     = 'broadband_cal';
        tpl.template_key       = sprintf('broadband_cal_%03d', j);
        tpl.band_list          = 1:B;
        tpl.band_id            = 0;                  % 0 表示全频带
        tpl.lambda_arrival     = Config.m0.source.lambda_broadband;
        tpl.life_mode          = BC.life_mode;
        tpl.life_param         = BC.life_param;

        % 功率：全频带统一功率（可扩展为每频带不同）
        tpl.tx_power_dBm       = BC.tx_power_dBm;
        tpl.tx_power_by_band   = BC.tx_power_dBm * ones(1, B);

        % 固定位置
        idx = randi(GridValid.Nvalid);
        tpl.fixed_pos_xy       = GridValid.xy(idx, :);

        SourceTemplates.broadband_cal(j) = tpl;
    end

    %% ===== 2. persistent_cal（持续存在的标校源，保底） =====
    PC = Config.m0.source.persistent_cal;
    tpl = struct();
    tpl.source_type     = 'persistent_cal';
    tpl.source_subtype  = 'persistent_cal';
    tpl.template_key    = 'persistent_cal_001';
    tpl.band_list       = PC.band_id;
    tpl.band_id         = PC.band_id;
    tpl.fixed_pos_xy    = PC.position;
    tpl.tx_power_dBm    = PC.tx_power_dBm;
    SourceTemplates.persistent_cal = tpl;

    %% ===== 3. target（待定位源，按频带） =====
    for b = 1:B
        Nocoop_b = Config.m0.target.Nocoop(b);
        for n = 1:Nocoop_b
            tpl = struct();
            tpl.source_type        = 'target';
            tpl.source_subtype     = 'target';
            tpl.template_key       = sprintf('target_b%d_%03d', b, n);
            tpl.band_list          = b;
            tpl.band_id            = b;
            tpl.lambda_arrival     = Config.m0.source.lambda_target(b);
            tpl.life_mode          = Config.m0.target.life_mode;
            tpl.life_param         = Config.m0.target.life_param;
            tpl.tx_power_range_dBm = Config.m0.target.power_range_dBm(b, :);
            tpl.position_mode      = Config.m0.target.position_mode;

            SourceTemplates.target(b, n) = tpl;
        end
    end

    %% ===== 4. opportunistic（机遇源，单一模板；实例运行时生成） =====
    tpl = struct();
    tpl.source_type    = 'opportunistic';
    tpl.source_subtype = 'opportunistic';
    tpl.template_key   = 'opportunistic_template';
    tpl.lambda_arrival = Config.m0.source.lambda_opportunistic;
    tpl.life_mode      = Config.m0.source.opportunistic.life_mode;
    tpl.life_param     = Config.m0.source.opportunistic.life_param;
    SourceTemplates.opportunistic = tpl;

    fprintf('[M0] 模板库: broadband_cal=%d, persistent_cal=1, target=%s, opportunistic=1\n', ...
        BC.count, mat2str(Config.m0.target.Nocoop));
end
