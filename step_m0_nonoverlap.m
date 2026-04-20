function [M0State, FrameState_t, M0Logs] = step_m0_nonoverlap( ...
    M0State, SourceTemplates, GridValid, Config, t, M0Logs)
% STEP_M0_NONOVERLAP  M0 单帧入口（统一四类源 + 优先级裁决）
%
%   流程：
%     1. update_persistent_cal_instances_m0  — always-on 保底
%     2. update_broadband_cal_instances_m0   — 宽带标校源（最高优先级）
%     3. update_opportunistic_instances_m0   — 机遇源（含三类位置先验）
%     4. update_target_instances_m0          — 待定位源（保留原逻辑）
%     5. resolve_band_occupancy_m0           — 优先级裁决（信源发生器内部）
%     6. finalize_frame_state_m0             — 组装统一 SourceEvent 输出
%     7. append_truth_logs_m0                — 真值日志

    M0State.frame_id = t;
    M0State.time_sec = t * Config.m0.dt;

    [M0State, PersistentCandidates] = update_persistent_cal_instances_m0( ...
        M0State, SourceTemplates, Config);

    [M0State, BroadbandCandidates] = update_broadband_cal_instances_m0( ...
        M0State, SourceTemplates, Config);

    [M0State, OppCandidatesPerBand] = update_opportunistic_instances_m0( ...
        M0State, SourceTemplates, GridValid, Config);

    [M0State, TargetCandidatesPerBand] = update_target_instances_m0( ...
        M0State, SourceTemplates, GridValid, Config);

    [M0State, BandWinners] = resolve_band_occupancy_m0( ...
        BroadbandCandidates, OppCandidatesPerBand, TargetCandidatesPerBand, ...
        PersistentCandidates, M0State, Config);

    FrameState_t = finalize_frame_state_m0(M0State, BandWinners, t, Config);

    M0Logs = append_truth_logs_m0(M0Logs, FrameState_t, Config);
end
