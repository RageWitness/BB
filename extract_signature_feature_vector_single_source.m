function [feat_vec, feat_names, match_weights] = ...
    extract_signature_feature_vector_single_source(tpl)
% EXTRACT_SIGNATURE_FEATURE_VECTOR_SINGLE_SOURCE  从模板构造静态特征向量
%
%   [feat_vec, feat_names, match_weights] =
%       extract_signature_feature_vector_single_source(tpl)
%
%   特征向量 12 维：
%     [band_mask(4), norm_power(1), stability(1), time_onehot(3), pos_onehot(3)]
%
%   输入：
%       tpl - 模板结构体（需含 band_mask, power_nominal_dBm,
%             power_stability_level, time_pattern_type, position_prior.mode）
%
%   输出：
%       feat_vec      - (1 x 12)
%       feat_names    - (1 x 12) cell
%       match_weights - (1 x 12) 匹配权重

    % 1. band_mask (4 维)
    bm = tpl.band_mask(:)';  % 确保 1 x B

    % 2. 归一化功率等级 (1 维)
    %    映射区间 [-20, 15] dBm → [0, 1]
    p_raw = tpl.power_nominal_dBm;
    p_norm = (p_raw - (-20)) / (15 - (-20));
    p_norm = max(0, min(1, p_norm));

    % 3. 功率稳定性等级 (1 维)
    %    映射 1/2/3 → 归一化到 [0, 1]
    s_norm = (tpl.power_stability_level - 1) / 2;

    % 4. 时间行为 one-hot (3 维)
    t_oh = encode_time_pattern_type(tpl.time_pattern_type);

    % 5. 位置先验 one-hot (3 维)
    p_oh = encode_position_prior_mode(tpl.position_prior.mode);

    % 组装
    feat_vec = [bm, p_norm, s_norm, t_oh, p_oh];

    feat_names = {'band1', 'band2', 'band3', 'band4', ...
                  'power_level', 'power_stability', ...
                  'time_scheduled', 'time_poisson', 'time_unknown', ...
                  'pos_fixed_known', 'pos_time_known', 'pos_unknown'};

    % 匹配权重：band_mask 和 position_prior 权重较高
    match_weights = [1, 1, 1, 1, ...    % band_mask
                     0.8, ...            % power_level
                     0.6, ...            % power_stability
                     0.7, 0.7, 0.7, ...  % time_pattern
                     1.0, 1.0, 1.0];     % position_prior
end
