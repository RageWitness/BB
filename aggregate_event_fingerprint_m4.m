function F_obs = aggregate_event_fingerprint_m4(event)
% AGGREGATE_EVENT_FINGERPRINT_M4  将事件窗口内观测聚合为单个 RSS 指纹
%
%   F_obs = aggregate_event_fingerprint_m4(event)
%
%   对事件 e 在频带 b 上的 obs_segment_dBm (M x L)，
%   取帧维度均值得到 M x 1 向量。
%
%   输入：
%       event - 单个事件结构体，需含 obs_segment_dBm (M x L)
%
%   输出：
%       F_obs - (M x 1) 事件平均 RSS 指纹 (dBm)

    seg = event.obs_segment_dBm;  % M x L

    if event.duration > 1
        F_obs = mean(seg, 2);  % M x 1
    else
        F_obs = seg(:);        % M x 1
    end
end
