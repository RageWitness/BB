function [F_obs_lin, F_obs_shape_l1] = aggregate_event_fingerprint_m4(event, agg_cfg)
% AGGREGATE_EVENT_FINGERPRINT_M4  将事件窗口内观测聚合为单个 RSS 指纹
%
%   [F_obs_lin, F_obs_shape_l1] = aggregate_event_fingerprint_m4(event)
%   [F_obs_lin, F_obs_shape_l1] = aggregate_event_fingerprint_m4(event, agg_cfg)
%
%   聚合模式：
%     'linear_mean'  — 线性域均值（默认）
%     'linear_median'— 线性域中位数
%
%   输出统一为线性域 (mW)，同时输出 L1 归一化形状向量。
%
%   输入：
%       event   - 单个事件结构体，需含 obs_segment_lin (M x L)
%       agg_cfg - [可选] 聚合配置
%                   agg_cfg.mode = 'linear_mean' / 'linear_median'
%
%   输出：
%       F_obs_lin       - (M x 1) 聚合后 RSS 指纹 (线性 mW)
%       F_obs_shape_l1  - (M x 1) L1 归一化形状向量

    % 默认配置
    if nargin < 2 || isempty(agg_cfg)
        agg_cfg = struct();
    end
    if ~isfield(agg_cfg, 'mode')
        agg_cfg.mode = 'linear_mean';
    end

    % 取线性域观测段
    if isfield(event, 'obs_segment_lin') && ~isempty(event.obs_segment_lin)
        seg_lin = event.obs_segment_lin;  % M x L
    else
        % 回退：从 dBm 转线性
        seg_lin = 10.^(event.obs_segment_dBm / 10);
    end

    % 单帧无需聚合
    if event.duration <= 1
        F_obs_lin = max(seg_lin(:), 1e-30);
        F_obs_shape_l1 = F_obs_lin / (sum(F_obs_lin) + 1e-30);
        return;
    end

    switch agg_cfg.mode
        case 'linear_mean'
            F_obs_lin = mean(seg_lin, 2);      % M x 1

        case 'linear_median'
            F_obs_lin = median(seg_lin, 2);    % M x 1

        otherwise
            error('[M4] 未知 agg_cfg.mode: %s', agg_cfg.mode);
    end

    % 确保非负
    F_obs_lin = max(F_obs_lin, 1e-30);

    % L1 归一化形状向量
    F_obs_shape_l1 = F_obs_lin / (sum(F_obs_lin) + 1e-30);
end
