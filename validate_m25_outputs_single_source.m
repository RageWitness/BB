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

        % L1 形状指纹维度检查
        if ~isfield(fb, 'F_shape_l1') || size(fb.F_shape_l1, 1) ~= M || size(fb.F_shape_l1, 2) ~= G
            fprintf('  [ERROR] Band %d: F_shape_l1 维度不正确或缺失\n', b);
            n_err = n_err + 1;
        end
        if ~isfield(fb, 'norm_l1') || length(fb.norm_l1) ~= G
            fprintf('  [ERROR] Band %d: norm_l1 长度不正确或缺失\n', b);
            n_err = n_err + 1;
        end
        % 形状指纹归一化检查：每列 L1 和应约为 1
        if isfield(fb, 'F_shape_l1')
            col_sums = sum(fb.F_shape_l1, 1);
            bad_cols = sum(abs(col_sums - 1) > 1e-6);
            if bad_cols > 0
                fprintf('  [WARN] Band %d: F_shape_l1 有 %d 列 L1 和偏离 1\n', b, bad_cols);
            end
        end

        % ---- 概率 shape 指纹校验 ----
        if isfield(fb, 'F_shape_prob_mu')
            if size(fb.F_shape_prob_mu, 1) ~= M || size(fb.F_shape_prob_mu, 2) ~= G
                fprintf('  [ERROR] Band %d: F_shape_prob_mu 维度 %dx%d, 应为 %dx%d\n', ...
                    b, size(fb.F_shape_prob_mu, 1), size(fb.F_shape_prob_mu, 2), M, G);
                n_err = n_err + 1;
            end
            if any(isnan(fb.F_shape_prob_mu(:))) || any(isinf(fb.F_shape_prob_mu(:)))
                fprintf('  [ERROR] Band %d: F_shape_prob_mu 含 NaN/Inf\n', b);
                n_err = n_err + 1;
            end
        end
        if isfield(fb, 'F_shape_prob_var')
            if size(fb.F_shape_prob_var, 1) ~= M || size(fb.F_shape_prob_var, 2) ~= G
                fprintf('  [ERROR] Band %d: F_shape_prob_var 维度 %dx%d, 应为 %dx%d\n', ...
                    b, size(fb.F_shape_prob_var, 1), size(fb.F_shape_prob_var, 2), M, G);
                n_err = n_err + 1;
            end
            if any(fb.F_shape_prob_var(:) < 0)
                fprintf('  [ERROR] Band %d: F_shape_prob_var 含负值\n', b);
                n_err = n_err + 1;
            end
            if any(isnan(fb.F_shape_prob_var(:))) || any(isinf(fb.F_shape_prob_var(:)))
                fprintf('  [ERROR] Band %d: F_shape_prob_var 含 NaN/Inf\n', b);
                n_err = n_err + 1;
            end
        end
        if isfield(fb, 'ap_weight_global')
            if length(fb.ap_weight_global) ~= M
                fprintf('  [ERROR] Band %d: ap_weight_global 长度 %d, 应为 %d\n', ...
                    b, length(fb.ap_weight_global), M);
                n_err = n_err + 1;
            end
            if any(isnan(fb.ap_weight_global)) || any(isinf(fb.ap_weight_global))
                fprintf('  [ERROR] Band %d: ap_weight_global 含 NaN/Inf\n', b);
                n_err = n_err + 1;
            end
            if any(fb.ap_weight_global <= 0)
                fprintf('  [WARN] Band %d: ap_weight_global 含非正值\n', b);
            end
        end

        % 指纹值范围检查（RSS 不应超出 [-200, 30] dBm）
        if isfield(fb, 'F_dBm')
            if any(fb.F_dBm(:) > 30) || any(fb.F_dBm(:) < -200)
                fprintf('  [WARN] Band %d: F_dBm 值范围 [%.1f, %.1f] 可能异常\n', ...
                    b, min(fb.F_dBm(:)), max(fb.F_dBm(:)));
            end
        end

        % ---- RF_minmax 校验（新框架主指纹） ----
        if isfield(fb, 'RF_minmax')
            if size(fb.RF_minmax, 1) ~= M || size(fb.RF_minmax, 2) ~= G
                fprintf('  [ERROR] Band %d: RF_minmax 维度 %dx%d, 应为 %dx%d\n', ...
                    b, size(fb.RF_minmax, 1), size(fb.RF_minmax, 2), M, G);
                n_err = n_err + 1;
            end
            if any(isnan(fb.RF_minmax(:))) || any(isinf(fb.RF_minmax(:)))
                fprintf('  [ERROR] Band %d: RF_minmax 含 NaN/Inf\n', b);
                n_err = n_err + 1;
            end
            mm_min = min(fb.RF_minmax(:));
            mm_max = max(fb.RF_minmax(:));
            if mm_min < -1.001 || mm_max > 1.001
                fprintf('  [WARN] Band %d: RF_minmax 范围 [%.4f, %.4f] 超出 [-1,1]\n', ...
                    b, mm_min, mm_max);
            end
            if isfield(fb, 'rf_minmax_n_degenerate') && fb.rf_minmax_n_degenerate > 0
                fprintf('  [INFO] Band %d: RF_minmax %d 个退化网格点（常值信号）\n', ...
                    b, fb.rf_minmax_n_degenerate);
            end
            fprintf('  [OK] Band %d: RF_minmax [%.4f, %.4f], mode=%s\n', ...
                b, mm_min, mm_max, fb.rf_norm_mode);
        end

        if isfield(fb, 'RF_raw')
            if size(fb.RF_raw, 1) ~= M || size(fb.RF_raw, 2) ~= G
                fprintf('  [ERROR] Band %d: RF_raw 维度不正确\n', b);
                n_err = n_err + 1;
            end
        end

        % ---- hard-negative 字段校验 ----
        if isfield(fb, 'hardneg_idx')
            if size(fb.hardneg_idx, 1) ~= G
                fprintf('  [ERROR] Band %d: hardneg_idx 行数 %d, 应为 %d\n', ...
                    b, size(fb.hardneg_idx, 1), G);
                n_err = n_err + 1;
            end
            if any(isnan(fb.hardneg_idx(:)))
                fprintf('  [ERROR] Band %d: hardneg_idx 含 NaN\n', b);
                n_err = n_err + 1;
            end
        end
        if isfield(fb, 'ap_weight_hn')
            if size(fb.ap_weight_hn, 1) ~= M || size(fb.ap_weight_hn, 2) ~= G
                fprintf('  [ERROR] Band %d: ap_weight_hn 维度 %dx%d, 应为 %dx%d\n', ...
                    b, size(fb.ap_weight_hn, 1), size(fb.ap_weight_hn, 2), M, G);
                n_err = n_err + 1;
            end
            if any(isnan(fb.ap_weight_hn(:))) || any(isinf(fb.ap_weight_hn(:)))
                fprintf('  [ERROR] Band %d: ap_weight_hn 含 NaN/Inf\n', b);
                n_err = n_err + 1;
            end
            if any(fb.ap_weight_hn(:) <= 0)
                fprintf('  [WARN] Band %d: ap_weight_hn 含非正值\n', b);
            end
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
