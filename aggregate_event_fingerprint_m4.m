function F_obs = aggregate_event_fingerprint_m4(event, agg_cfg)
% AGGREGATE_EVENT_FINGERPRINT_M4  将事件窗口内观测聚合为单个 RSS 指纹
%
%   F_obs = aggregate_event_fingerprint_m4(event)
%   F_obs = aggregate_event_fingerprint_m4(event, agg_cfg)
%
%   三种聚合模式：
%     'linear_mean' — 线性域均值后转 dBm（默认，对对数正态更无偏）
%     'dbm_median'  — dBm 域逐 AP 中位数（鲁棒）
%     'dbm_mean'    — dBm 域均值（旧行为）
%
%   输入：
%       event   - 单个事件结构体
%       agg_cfg - [可选] 聚合配置
%                   agg_cfg.mode = 'linear_mean' / 'dbm_median' / 'dbm_mean'
%
%   输出：
%       F_obs   - (M x 1) 聚合后 RSS 指纹 (dBm)

    % 默认配置
    if nargin < 2 || isempty(agg_cfg)
        agg_cfg = struct();
    end
    if ~isfield(agg_cfg, 'mode')
        agg_cfg.mode = 'linear_mean';
    end

    % 单帧事件无需聚合
    if event.duration <= 1
        if isfield(event, 'obs_segment_dBm')
            F_obs = event.obs_segment_dBm(:);
        else
            F_obs = 10 * log10(event.obs_segment_lin(:));
        end
        return;
    end

    switch agg_cfg.mode
        case 'linear_mean'
            if isfield(event, 'obs_segment_lin') && ~isempty(event.obs_segment_lin)
                seg_lin = event.obs_segment_lin;           % M x L
                F_lin = mean(seg_lin, 2);                  % M x 1
                F_obs = 10 * log10(max(F_lin, 1e-20));
            else
                % 回退：从 dBm 转线性再取均值
                seg_dBm = event.obs_segment_dBm;           % M x L
                seg_lin = 10.^(seg_dBm / 10);
                F_lin = mean(seg_lin, 2);
                F_obs = 10 * log10(max(F_lin, 1e-20));
            end

        case 'dbm_median'
            seg_dBm = event.obs_segment_dBm;               % M x L
            F_obs = median(seg_dBm, 2);                    % M x 1

        case 'dbm_mean'
            seg_dBm = event.obs_segment_dBm;               % M x L
            F_obs = mean(seg_dBm, 2);                      % M x 1

        otherwise
            error('[M4] 未知 agg_cfg.mode: %s', agg_cfg.mode);
    end
end
