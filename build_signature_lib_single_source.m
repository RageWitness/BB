function SignatureLib = build_signature_lib_single_source(SourceTemplates, Config)
% BUILD_SIGNATURE_LIB_SINGLE_SOURCE  从 SourceTemplates 构建静态模板库
%
%   SignatureLib = build_signature_lib_single_source(SourceTemplates, Config)
%
%   将 SourceTemplates 中四类模板展平为统一数组，
%   提取静态特征字段并预计算特征向量。
%
%   注意：只存静态模板特征，不存在线可信度证据。
%
%   输出：
%       SignatureLib.n_templates
%       SignatureLib.templates(j) — 统一索引的模板特征数组

    B = Config.m0.num_bands;
    idx = 0;  % 模板计数器

    %% 1. trusted_fixed
    for j = 1:numel(SourceTemplates.trusted)
        tpl = SourceTemplates.trusted(j);
        idx = idx + 1;
        templates(idx) = extract_template_entry(tpl);
    end

    %% 2. prior_pos_known
    for b = 1:B
        Sprior_b = Config.m0.prior.Sprior(b);
        for s = 1:Sprior_b
            tpl = SourceTemplates.prior_pos(b, s);
            idx = idx + 1;
            templates(idx) = extract_template_entry(tpl);
        end
    end

    %% 3. prior_time_known
    for b = 1:B
        Tprior_b = Config.m0.prior.Tprior(b);
        for q = 1:Tprior_b
            tpl = SourceTemplates.prior_time(b, q);
            idx = idx + 1;
            templates(idx) = extract_template_entry(tpl);
        end
    end

    %% 4. ordinary_target
    for b = 1:B
        Nocoop_b = Config.m0.target.Nocoop(b);
        for n = 1:Nocoop_b
            tpl = SourceTemplates.target(b, n);
            idx = idx + 1;
            templates(idx) = extract_template_entry(tpl);
        end
    end

    SignatureLib.n_templates = idx;
    SignatureLib.templates   = templates;

    fprintf('[M2.5] SignatureLib 构建完成: %d 个模板\n', idx);

    % 打印摘要
    types = {templates.source_type};
    subtypes = {templates.source_subtype};
    n_trusted = sum(strcmp(types, 'trusted_fixed'));
    n_pp = sum(strcmp(subtypes, 'prior_pos_known'));
    n_pt = sum(strcmp(subtypes, 'prior_time_known'));
    n_tg = sum(strcmp(types, 'ordinary_target'));
    fprintf('  trusted_fixed=%d, prior_pos_known=%d, prior_time_known=%d, ordinary_target=%d\n', ...
        n_trusted, n_pp, n_pt, n_tg);
end


function entry = extract_template_entry(tpl)
% EXTRACT_TEMPLATE_ENTRY  从一个模板中提取 SignatureLib 条目

    % 基本标识
    entry.source_type        = tpl.source_type;
    entry.source_subtype     = tpl.source_subtype;
    entry.template_namespace = tpl.template_namespace;
    entry.template_key       = tpl.template_key;

    % 频带覆盖
    entry.band_mask          = tpl.band_mask;

    % 功率特征
    entry.power_nominal_dBm    = tpl.power_nominal_dBm;
    entry.power_range_dBm      = tpl.power_range_dBm;
    entry.power_stability_level = tpl.power_stability_level;

    % 时间行为
    entry.time_pattern_type     = tpl.time_pattern_type;
    entry.expected_duration_range = tpl.expected_duration_range;

    % 位置先验（对象化，不用布尔）
    entry.known_position_mode  = tpl.position_prior.mode;
    if isfield(tpl.position_prior, 'candidate_positions') && ...
            ~isempty(tpl.position_prior.candidate_positions)
        entry.candidate_positions = tpl.position_prior.candidate_positions;
    else
        entry.candidate_positions = [];
    end

    % 可信度
    entry.credibility_prior    = tpl.credibility_prior;

    % 预计算特征向量
    [entry.feat_vec, entry.feat_names, entry.match_weights] = ...
        extract_signature_feature_vector_single_source(tpl);
end
