function [D_vec, extra] = compute_wknn_distance_fp_l2( ...
    F_obs_lin, F_obs_shape_l1, F_obs_dBm, SpatialFP_band, match_cfg)
% COMPUTE_WKNN_DISTANCE_FP_L2  指纹自身域上的 L1/L2 WKNN 距离
%
%   支持的 fingerprint_type:
%     'rf_minmax'    - Min-Max 标准化到 [-1,1] (库 + 观测同步标准化)
%     'rf_raw'       - 原始线性功率
%     'shape_l1'     - L1 归一化的 shape 向量
%     'centered_dBm' - dBm 域去均值后的中心化向量
%
%   match_cfg.fp_distance: 'L1' 或 'L2'

    fp_type = match_cfg.fingerprint_type;
    fp_dist = upper(match_cfg.fp_distance);

    switch fp_type
        case 'rf_minmax'
            F_lib = SpatialFP_band.RF_minmax;          % M x G_sub
            F_obs = minmax_normalize_obs(F_obs_lin);   % M x 1

        case 'rf_raw'
            if isfield(SpatialFP_band, 'RF_raw')
                F_lib = SpatialFP_band.RF_raw;
            else
                F_lib = SpatialFP_band.F_lin;
            end
            F_obs = F_obs_lin(:);

        case 'shape_l1'
            F_lib = SpatialFP_band.F_shape_l1;
            F_obs = F_obs_shape_l1(:);

        case 'centered_dBm'
            F_lib = SpatialFP_band.centered_dBm;
            F_obs = F_obs_dBm(:) - mean(F_obs_dBm(:));

        otherwise
            error('[M4] compute_wknn_distance_fp_l2: 未知 fingerprint_type=%s', fp_type);
    end

    if size(F_obs, 1) ~= size(F_lib, 1)
        error('[M4] compute_wknn_distance_fp_l2: F_obs (M=%d) 与库 (M=%d) AP 数不匹配', ...
            size(F_obs, 1), size(F_lib, 1));
    end

    diff = F_obs - F_lib;   % 利用广播 M x G_sub
    switch fp_dist
        case 'L2'
            D_vec = sqrt(sum(diff.^2, 1));
        case 'L1'
            D_vec = sum(abs(diff), 1);
        otherwise
            error('[M4] compute_wknn_distance_fp_l2: 未知 fp_distance=%s', fp_dist);
    end

    extra = struct();
    extra.fingerprint_type = fp_type;
    extra.fp_distance      = fp_dist;
    extra.F_obs_used       = F_obs;
end


function v = minmax_normalize_obs(F_obs_lin)
    F_obs_lin = F_obs_lin(:);
    z_min = min(F_obs_lin);
    z_max = max(F_obs_lin);
    z_range = z_max - z_min;
    if z_range < 1e-30
        v = zeros(size(F_obs_lin));
    else
        v = 2 * (F_obs_lin - z_min) / z_range - 1;
    end
end
