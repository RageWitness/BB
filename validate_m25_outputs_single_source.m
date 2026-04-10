function validate_m25_outputs_single_source(SpatialFP, SignatureLib)
% VALIDATE_M25_OUTPUTS_SINGLE_SOURCE  校验 M2.5 输出的完整性和一致性
%
%   validate_m25_outputs_single_source(SpatialFP, SignatureLib)
%
%   检查项：
%     1. SpatialFP 维度合法
%     2. SignatureLib 字段齐全
%     3. template_key 唯一
%     4. 指纹值范围合理

    fprintf('[M2.5] 开始校验...\n');
    n_err = 0;

    %% 1. SpatialFP 维度
    B = SpatialFP.B;
    M = SpatialFP.M;
    G = SpatialFP.G;

    assert(size(SpatialFP.grid_xy, 1) == G && size(SpatialFP.grid_xy, 2) == 2, ...
        'grid_xy 维度不正确');

    for b = 1:B
        fb = SpatialFP.band(b);
        if size(fb.F_dBm, 1) ~= M || size(fb.F_dBm, 2) ~= G
            fprintf('  [ERROR] Band %d: F_dBm 维度 %dx%d, 应为 %dx%d\n', ...
                b, size(fb.F_dBm, 1), size(fb.F_dBm, 2), M, G);
            n_err = n_err + 1;
        end
        if size(fb.centered_dBm, 1) ~= M || size(fb.centered_dBm, 2) ~= G
            fprintf('  [ERROR] Band %d: centered_dBm 维度不正确\n', b);
            n_err = n_err + 1;
        end
        if length(fb.mean_dBm) ~= G
            fprintf('  [ERROR] Band %d: mean_dBm 长度不正确\n', b);
            n_err = n_err + 1;
        end

        % 指纹值范围检查（RSS 不应超出 [-200, 30] dBm）
        if any(fb.F_dBm(:) > 30) || any(fb.F_dBm(:) < -200)
            fprintf('  [WARN] Band %d: F_dBm 值范围 [%.1f, %.1f] 可能异常\n', ...
                b, min(fb.F_dBm(:)), max(fb.F_dBm(:)));
        end
    end

    %% 2. SignatureLib 字段检查
    required_fields = {'source_type', 'source_subtype', 'template_namespace', ...
        'template_key', 'band_mask', 'power_nominal_dBm', 'power_range_dBm', ...
        'power_stability_level', 'time_pattern_type', 'expected_duration_range', ...
        'known_position_mode', 'candidate_positions', 'credibility_prior', ...
        'feat_vec', 'feat_names', 'match_weights'};

    for fi = 1:numel(required_fields)
        fn = required_fields{fi};
        if ~isfield(SignatureLib.templates, fn)
            fprintf('  [ERROR] SignatureLib 缺少字段: %s\n', fn);
            n_err = n_err + 1;
        end
    end

    %% 3. template_key 唯一性
    keys = {SignatureLib.templates.template_key};
    if numel(unique(keys)) ~= numel(keys)
        fprintf('  [ERROR] template_key 存在重复\n');
        n_err = n_err + 1;
    end

    %% 4. feat_vec 维度一致性
    for j = 1:SignatureLib.n_templates
        if length(SignatureLib.templates(j).feat_vec) ~= 12
            fprintf('  [ERROR] 模板 %s 的 feat_vec 长度 %d, 应为 12\n', ...
                SignatureLib.templates(j).template_key, ...
                length(SignatureLib.templates(j).feat_vec));
            n_err = n_err + 1;
        end
    end

    if n_err == 0
        fprintf('[M2.5] [PASS] 全部校验通过\n');
    else
        fprintf('[M2.5] [FAIL] 发现 %d 个错误\n', n_err);
    end
end
