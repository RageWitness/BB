function [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
    M0State, SourceTemplates, GridValid, Config, t, M0Logs)
% STEP_M0_NONOVERLAP  M0 单帧更新总入口
%
%   [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap(...)
%
%   内部流程：
%     1. update_trusted_instances_m0  — 更新 trusted_fixed 到达与寿命
%     2. update_prior_instances_m0    — 更新 prior 两类源的 schedule
%     3. update_target_instances_m0   — 更新 ordinary_target 到达与寿命
%     4. resolve_band_occupancy_m0    — 裁决每频带唯一主导源
%     5. finalize_frame_state_m0      — 组装帧输出
%     6. append_truth_logs_m0         — 记录真值日志
%
%   输入：
%       M0State         - 运行时状态
%       SourceTemplates - 模板库
%       GridValid       - 有效网格
%       Config          - 配置
%       t               - 当前帧号 (1-indexed)
%       M0Logs          - 真值日志
%
%   输出：
%       M0State         - 更新后的运行时状态
%       FrameState_t    - 当前帧输出
%       M0Logs          - 更新后的真值日志

    % 更新帧号与时间
    M0State.frame_id = t;
    M0State.time_sec = t * Config.m0.dt;

    %% 1. 更新 trusted_fixed 实例
    [M0State, TrustedCandidates] = update_trusted_instances_m0( ...
        M0State, SourceTemplates, Config, t);

    %% 2. 更新 prior 实例
    [M0State, PriorCandidatesPerBand] = update_prior_instances_m0( ...
        M0State, SourceTemplates, GridValid, Config, t);

    %% 3. 更新 ordinary_target 实例
    [M0State, TargetCandidatesPerBand] = update_target_instances_m0( ...
        M0State, SourceTemplates, GridValid, Config, t);

    %% 4. 频带占用裁决
    [M0State, BandWinners] = resolve_band_occupancy_m0( ...
        TrustedCandidates, PriorCandidatesPerBand, TargetCandidatesPerBand, ...
        M0State, Config);

    %% 5. 组装帧输出
    FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config);

    %% 6. 记录真值日志
    M0Logs = append_truth_logs_m0(M0Logs, FrameState_t, Config);
end
