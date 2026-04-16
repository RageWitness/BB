function GroupList = group_events_multiband_m3(EventList, Config)
% GROUP_EVENTS_MULTIBAND_M3  跨频带事件关联（按时间重叠/间隔建组）
%
%   GroupList = group_events_multiband_m3(EventList, Config)
%
%   规则：
%     1) 时间重叠或间隔 <= max_gap_frames 的事件可进入同一 group
%     2) group 内每个 band 只保留一个成员事件
%     3) 同 band 冲突时，优先保留与 group 时间跨度更一致的事件

    GroupList = [];
    if isempty(EventList)
        return;
    end

    cfg = fill_group_defaults(Config);
    B = max(cfg.num_bands, max([EventList.band_id]));
    n_events = numel(EventList);

    % 按起始时间排序，逐个并入已有 group
    [~, order] = sortrows([[EventList.t_start].', [EventList.t_end].', [EventList.band_id].'], [1, 2, 3]); %#ok<NBRAK2>
    sorted_idx = order(:).';

    groups = [];
    for ii = 1:numel(sorted_idx)
        eidx = sorted_idx(ii);
        ev = EventList(eidx);

        [cand_gid, cand_cost] = find_candidate_groups(ev, groups, cfg.max_gap_frames);
        if isempty(cand_gid)
            groups = append_new_group(groups, ev, eidx, B, EventList);
            continue;
        end

        [~, rank_idx] = sort(cand_cost, 'ascend');
        assigned = false;
        for kk = 1:numel(rank_idx)
            gid = cand_gid(rank_idx(kk));
            b = ev.band_id;
            old_idx = groups(gid).member_event_idx_per_band(b);
            if old_idx == 0
                groups(gid).member_event_idx_per_band(b) = eidx;
                groups(gid) = refresh_group_stats(groups(gid), EventList, B);
                assigned = true;
                break;
            end

            old_ev = EventList(old_idx);
            old_score = temporal_consistency_score(old_ev, groups(gid));
            new_score = temporal_consistency_score(ev, groups(gid));
            if new_score > old_score + cfg.same_band_replace_margin
                groups(gid).member_event_idx_per_band(b) = eidx;
                groups(gid) = refresh_group_stats(groups(gid), EventList, B);
                assigned = true;
                break;
            end
        end

        if ~assigned
            groups = append_new_group(groups, ev, eidx, B, EventList);
        end
    end

    % 重排 group_id（按时间）
    if isempty(groups)
        GroupList = [];
        return;
    end
    [~, gorder] = sortrows([[groups.start_frame].', [groups.end_frame].'], [1, 2]); %#ok<NBRAK2>
    groups = groups(gorder);
    for g = 1:numel(groups)
        groups(g).group_id = g;
    end

    GroupList = groups;
    fprintf('[M3] group 构建完成: %d groups from %d events\n', numel(GroupList), n_events);
end


%% ==================== 局部函数 ====================

function cfg = fill_group_defaults(Config)
    cfg.max_gap_frames = 2;
    cfg.same_band_replace_margin = 0.05;
    cfg.num_bands = 4;

    if isfield(Config, 'm0') && isfield(Config.m0, 'num_bands')
        cfg.num_bands = Config.m0.num_bands;
    end
    if isfield(Config, 'm3') && isfield(Config.m3, 'grouping')
        g = Config.m3.grouping;
        if isfield(g, 'max_gap_frames')
            cfg.max_gap_frames = g.max_gap_frames;
        end
        if isfield(g, 'same_band_replace_margin')
            cfg.same_band_replace_margin = g.same_band_replace_margin;
        end
    end
end


function [cand_gid, cand_cost] = find_candidate_groups(ev, groups, max_gap_frames)
    cand_gid = [];
    cand_cost = [];
    for g = 1:numel(groups)
        gap = interval_gap_frames(ev.t_start, ev.t_end, groups(g).start_frame, groups(g).end_frame);
        if gap <= max_gap_frames
            c = temporal_attach_cost(ev, groups(g));
            cand_gid(end+1) = g; %#ok<AGROW>
            cand_cost(end+1) = c; %#ok<AGROW>
        end
    end
end


function gap = interval_gap_frames(s1, e1, s2, e2)
    if e1 < s2
        gap = s2 - e1 - 1;
    elseif e2 < s1
        gap = s1 - e2 - 1;
    else
        gap = 0;
    end
end


function c = temporal_attach_cost(ev, grp)
% 越小越好：优先时间跨度扩展小、中心偏差小
    new_start = min(ev.t_start, grp.start_frame);
    new_end = max(ev.t_end, grp.end_frame);
    old_span = grp.end_frame - grp.start_frame + 1;
    new_span = new_end - new_start + 1;
    span_expand = new_span - old_span;

    ev_center = 0.5 * (ev.t_start + ev.t_end);
    grp_center = 0.5 * (grp.start_frame + grp.end_frame);
    center_diff = abs(ev_center - grp_center);

    c = span_expand + 0.25 * center_diff;
end


function s = temporal_consistency_score(ev, grp)
% 越大越好
    ov = overlap_ratio(ev.t_start, ev.t_end, grp.start_frame, grp.end_frame);
    len_ratio = min(ev.duration, grp.duration_frames) / max(ev.duration, grp.duration_frames);
    s = 0.7 * ov + 0.3 * len_ratio;
end


function r = overlap_ratio(s1, e1, s2, e2)
    ov_start = max(s1, s2);
    ov_end = min(e1, e2);
    if ov_end < ov_start
        r = 0;
        return;
    end
    ov_len = ov_end - ov_start + 1;
    union_len = max(e1, e2) - min(s1, s2) + 1;
    r = ov_len / max(union_len, 1);
end


function groups = append_new_group(groups, ev, eidx, B, EventList)
    grp = struct();
    grp.group_id = numel(groups) + 1;
    grp.start_frame = ev.t_start;
    grp.end_frame = ev.t_end;
    grp.duration_frames = ev.duration;
    grp.band_mask = zeros(1, B);
    grp.active_bands = ev.band_id;
    grp.active_band_count = 1;

    grp.member_event_idx_per_band = zeros(1, B);
    grp.member_event_idx_per_band(ev.band_id) = eidx;
    grp.member_event_indices = eidx;
    grp.member_event_ids = ev.event_id;

    if isempty(groups)
        groups = grp;
    else
        groups(end+1) = grp; %#ok<AGROW>
    end
    groups(end) = refresh_group_stats(groups(end), EventList, B);
end


function grp = refresh_group_stats(grp, EventList, B)
    idx = grp.member_event_idx_per_band;
    idx = idx(idx > 0);

    if isempty(idx)
        grp.start_frame = 0;
        grp.end_frame = 0;
        grp.duration_frames = 0;
        grp.band_mask = zeros(1, B);
        grp.active_bands = [];
        grp.active_band_count = 0;
        grp.member_event_indices = [];
        grp.member_event_ids = [];
        return;
    end

    starts = [EventList(idx).t_start];
    ends = [EventList(idx).t_end];
    bands = [EventList(idx).band_id];
    ids = [EventList(idx).event_id];

    grp.start_frame = min(starts);
    grp.end_frame = max(ends);
    grp.duration_frames = grp.end_frame - grp.start_frame + 1;

    grp.band_mask = zeros(1, B);
    grp.band_mask(bands) = 1;
    grp.active_bands = find(grp.band_mask > 0);
    grp.active_band_count = numel(grp.active_bands);
    grp.member_event_indices = idx;
    grp.member_event_ids = ids;
end
