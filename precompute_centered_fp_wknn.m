function SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band, fp_mode)
% PRECOMPUTE_CENTERED_FP_WKNN  预计算 WKNN 所需的指纹表示
%
%   SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band)
%   SpatialFP_band = precompute_centered_fp_wknn(SpatialFP_band, fp_mode)
%
%   fp_mode:
%     'rf_minmax'  — 新框架主模式：Min-Max 标准化到 [-1, 1]
%     'legacy'     — 旧模式：dBm centered + L1 shape（过渡兼容）
%     'both'       — 同时计算（默认）
%
%   === 新框架主指纹 ===
%   RF_raw(:,g)     = F_lin(:,g)                        原始线性功率
%   RF_minmax(:,g)  = 2*(z - Z_min)/(Z_max - Z_min) - 1 标准化到 [-1,1]
%   rf_norm_mode    = 'minmax_neg1_pos1'
%
%   === Legacy 指纹（过渡兼容，不再是主入口）===
%   mean_dBm, centered_dBm, mean_lin, centered_lin, norm_l1, F_shape_l1

    if nargin < 2 || isempty(fp_mode)
        fp_mode = 'both';
    end

    F_lin = SpatialFP_band.F_lin;   % M x G

    % ============================================
    %  新框架主模式：Min-Max 标准化到 [-1, 1]
    % ============================================
    if strcmp(fp_mode, 'rf_minmax') || strcmp(fp_mode, 'both')
        SpatialFP_band.RF_raw = SpatialFP_band.F_dBm;

        z_min = min(F_lin, [], 1);   % 1 x G
        z_max = max(F_lin, [], 1);   % 1 x G
        z_range = z_max - z_min;     % 1 x G

        % 退化处理：Z_max == Z_min 时输出全 0
        degenerate = z_range < 1e-30;
        z_range_safe = z_range;
        z_range_safe(degenerate) = 1;  % 防除零

        RF_minmax = 2 * (F_lin - z_min) ./ z_range_safe - 1;
        RF_minmax(:, degenerate) = 0;

        SpatialFP_band.RF_minmax = RF_minmax;
        SpatialFP_band.rf_norm_mode = 'minmax_neg1_pos1';
        SpatialFP_band.rf_minmax_zmin = z_min;
        SpatialFP_band.rf_minmax_zmax = z_max;
        SpatialFP_band.rf_minmax_n_degenerate = sum(degenerate);

        if any(degenerate)
            fprintf('    [MinMax] %d / %d 网格点为常值（退化为 0）\n', sum(degenerate), numel(degenerate));
        end
    end

    % ============================================
    %  Legacy 模式（过渡兼容，不再是主指纹构建入口）
    % ============================================
    if strcmp(fp_mode, 'legacy') || strcmp(fp_mode, 'both')
        % --- dBm 域（legacy） ---
        F_dBm = SpatialFP_band.F_dBm;      % M x G
        mu_dBm = mean(F_dBm, 1);           % 1 x G
        SpatialFP_band.mean_dBm     = mu_dBm;
        SpatialFP_band.centered_dBm = F_dBm - mu_dBm;

        % --- 线性域（legacy） ---
        mu_lin = mean(F_lin, 1);           % 1 x G
        mu_lin_safe = max(mu_lin, 1e-30);
        SpatialFP_band.mean_lin     = mu_lin;
        SpatialFP_band.centered_lin = F_lin ./ mu_lin_safe;

        % --- L1 shape（legacy） ---
        norm_l1 = sum(F_lin, 1);
        norm_l1_safe = max(norm_l1, 1e-30);
        SpatialFP_band.norm_l1    = norm_l1;
        SpatialFP_band.F_shape_l1 = F_lin ./ norm_l1_safe;

        SpatialFP_band.shape_library_mode = 'deterministic_l1_from_mean_lin';
    end
end
