function SignatureLib = build_signature_lib_single_source(SourceTemplates, Config)
% BUILD_SIGNATURE_LIB_SINGLE_SOURCE  从 SourceTemplates 构建静态模板库（适配新四类源）
%
%   将 broadband_cal / persistent_cal / target / opportunistic 四类模板
%   展平为统一数组，提取静态特征字段并预计算特征向量。

    B = Config.m0.num_bands;
    idx = 0;

    %% 1. broadband_cal
    if isfield(SourceTemplates, 'broadband_cal')
        for j = 1:numel(SourceTemplates.broadband_cal)
            idx = idx + 1;
            templates(idx) = extract_template_entry( ...
                SourceTemplates.broadband_cal(j), B);
        end
    end

    %% 2. persistent_cal
    if isfield(SourceTemplates, 'persistent_cal')
        idx = idx + 1;
        templates(idx) = extract_template_entry(SourceTemplates.persistent_cal, B);
    end

    %% 3. target
    if isfield(SourceTemplates, 'target')
        for b = 1:B
            Nocoop_b = Config.m0.target.Nocoop(b);
            for n = 1:Nocoop_b
                idx = idx + 1;
                templates(idx) = extract_template_entry( ...
                    SourceTemplates.target(b, n), B);
            end
        end
    end

    %% 4. opportunistic
    if isfield(SourceTemplates, 'opportunistic')
        idx = idx + 1;
        templates(idx) = extract_template_entry(SourceTemplates.opportunistic, B);
    end

    SignatureLib.n_templates = idx;
    SignatureLib.templates   = templates;

    fprintf('[M2.5] SignatureLib: %d 个模板 (broadband_cal/persistent_cal/target/opportunistic)\n', idx);
end


function entry = extract_template_entry(tpl, B)
    entry.source_type        = getf(tpl, 'source_type', '');
    entry.source_subtype     = getf(tpl, 'source_subtype', entry.source_type);
    entry.template_namespace = getf(tpl, 'template_namespace', entry.source_type);
    entry.template_key       = getf(tpl, 'template_key', '');

    % 频带覆盖
    if isfield(tpl, 'band_mask') && ~isempty(tpl.band_mask)
        entry.band_mask = tpl.band_mask;
    elseif isfield(tpl, 'band_list') && ~isempty(tpl.band_list)
        entry.band_mask = zeros(1, B);
        valid = tpl.band_list(tpl.band_list >= 1 & tpl.band_list <= B);
        entry.band_mask(valid) = 1;
    else
        entry.band_mask = zeros(1, B);
    end

    % 功率特征
    if isfield(tpl, 'tx_power_dBm') && ~isempty(tpl.tx_power_dBm)
        entry.power_nominal_dBm = tpl.tx_power_dBm;
        entry.power_range_dBm   = [tpl.tx_power_dBm, tpl.tx_power_dBm];
    elseif isfield(tpl, 'tx_power_range_dBm')
        entry.power_nominal_dBm = mean(tpl.tx_power_range_dBm);
        entry.power_range_dBm   = tpl.tx_power_range_dBm;
    else
        entry.power_nominal_dBm = 0;
        entry.power_range_dBm   = [0 0];
    end
    entry.power_stability_level = stability_for(entry.source_type);

    % 时间行为
    entry.time_pattern_type       = 'poisson';
    entry.expected_duration_range = [5, 30];

    % 位置先验
    if isfield(tpl, 'fixed_pos_xy') && ~isempty(tpl.fixed_pos_xy)
        entry.known_position_mode  = 'fixed_known';
        entry.candidate_positions  = tpl.fixed_pos_xy;
    else
        entry.known_position_mode  = 'unknown';
        entry.candidate_positions  = [];
    end

    entry.credibility_prior = credibility_for(entry.source_type);

    [entry.feat_vec, entry.feat_names, entry.match_weights] = ...
        extract_signature_feature_vector_single_source(template_with_priors(tpl, entry));
end


function v = getf(s, fn, default)
    if isfield(s, fn) && ~isempty(s.(fn))
        v = s.(fn);
    else
        v = default;
    end
end


function lvl = stability_for(src_type)
    switch src_type
        case {'broadband_cal','persistent_cal'}, lvl = 3;
        case 'opportunistic',                    lvl = 2;
        otherwise,                                lvl = 1;
    end
end


function c = credibility_for(src_type)
    switch src_type
        case {'broadband_cal','persistent_cal'}, c = 0.99;
        case 'opportunistic',                    c = 0.75;
        otherwise,                                c = 0.2;
    end
end


function tpl = template_with_priors(tpl, entry)
% 给老的 extract_signature_feature_vector_single_source 提供它需要的字段
    if ~isfield(tpl, 'position_prior')
        tpl.position_prior = struct('mode', entry.known_position_mode, ...
            'candidate_positions', entry.candidate_positions, 'match_radius', []);
    end
    if ~isfield(tpl, 'time_prior')
        tpl.time_prior = struct('mode','unknown','schedule',[], ...
            'period_frames',[],'duration_frames',[],'phase_frames',[]);
    end
    if ~isfield(tpl, 'power_nominal_dBm'), tpl.power_nominal_dBm = entry.power_nominal_dBm; end
    if ~isfield(tpl, 'power_range_dBm'),   tpl.power_range_dBm = entry.power_range_dBm;     end
    if ~isfield(tpl, 'power_stability_level'), tpl.power_stability_level = entry.power_stability_level; end
    if ~isfield(tpl, 'credibility_prior'), tpl.credibility_prior = entry.credibility_prior; end
    if ~isfield(tpl, 'upgrade_potential'), tpl.upgrade_potential = 'none'; end
    if ~isfield(tpl, 'time_pattern_type'), tpl.time_pattern_type = entry.time_pattern_type; end
    if ~isfield(tpl, 'expected_duration_range'), tpl.expected_duration_range = entry.expected_duration_range; end
    if ~isfield(tpl, 'band_mask'), tpl.band_mask = entry.band_mask; end
end
