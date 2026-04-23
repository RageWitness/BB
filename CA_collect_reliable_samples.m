function CASamples = CA_collect_reliable_samples(EventList, Config)
% CA_COLLECT_RELIABLE_SAMPLES  从 EventList 筛选 CA 可靠样本
%
%   只保留 label==1 (persistent_cal) 或 label==2 (broadband_cal)，
%   且满足 exact location + exact power 的样本。

    Config = CA_fill_defaults(Config);
    LC = source_label_constants();

    CASamples = [];
    sid = 0;

    for i = 1:numel(EventList)
        ev = EventList(i);

        if ev.label ~= LC.PERSISTENT_CAL && ev.label ~= LC.BROADBAND_CAL
            continue;
        end

        sid = sid + 1;
        s = struct();
        s.sample_id  = sid;
        s.event_id   = ev.event_id;
        s.label      = ev.label;
        s.type_hat   = ev.type_hat;
        s.band_id    = ev.band_id;
        s.time_range = ev.time_range;
        s.duration   = ev.duration;
        s.source_uid = ev.source_uid;
        s.linked_template_key = ev.linked_template_key;
        s.location_prior = ev.location_prior;
        s.power_prior    = ev.power_prior;

        if isfield(ev, 'band_coverage_vec')
            s.band_list = find(ev.band_coverage_vec > 0);
        else
            s.band_list = ev.band_id;
        end
        if isfield(ev, 'instance_id')
            s.instance_id = ev.instance_id;
        else
            s.instance_id = 0;
        end

        s.obs_segment_lin = ev.obs_segment_lin;
        s.obs_segment_dBm = ev.obs_segment_dBm;

        s.valid = true;
        s.reject_reason = '';

        % --- 维度检查 ---
        [M_obs, ~] = size(s.obs_segment_lin);
        if M_obs < 2
            s.valid = false;
            s.reject_reason = 'obs_dimension_invalid';
        end

        % --- location 检查 ---
        if s.valid
            if ~isfield(s.location_prior, 'type') || ...
               ~strcmp(s.location_prior.type, 'exact')
                s.valid = false;
                s.reject_reason = 'location_not_exact';
            end
        end

        % --- power 检查 ---
        if s.valid
            pt = s.power_prior.type;
            if strcmp(pt, 'exact')
                % ok
            elseif strcmp(pt, 'exact_by_band')
                if ev.band_id > numel(s.power_prior.value)
                    s.valid = false;
                    s.reject_reason = 'power_not_exact_for_band';
                end
            else
                s.valid = false;
                s.reject_reason = 'power_not_exact';
            end
        end

        % --- position ---
        if s.valid && ~isempty(s.location_prior.value)
            s.position_xy = s.location_prior.value(:)';
        else
            s.position_xy = [];
            if s.valid
                s.valid = false;
                s.reject_reason = 'location_value_empty';
            end
        end

        % --- centered observation ---
        s.obs_centered_dBm = [];
        s.obs_mean_dBm     = [];
        s.obs_mean_lin     = [];
        s.valid_ap_mask    = [];

        if s.valid
            [z_c, y_dBm, y_lin, vam, cinfo] = ...
                CA_compute_centered_observation(s.obs_segment_lin, Config);

            s.obs_centered_dBm = z_c;
            s.obs_mean_dBm     = y_dBm;
            s.obs_mean_lin     = y_lin;
            s.valid_ap_mask    = vam;

            if strcmp(cinfo.status, 'too_few_valid_aps')
                s.valid = false;
                s.reject_reason = 'too_few_valid_aps';
            end
        end

        if isempty(CASamples)
            CASamples = s;
        else
            CASamples(end+1) = s; %#ok<AGROW>
        end
    end

    if isempty(CASamples)
        CASamples = struct([]);
    end

    n_total = numel(CASamples);
    if n_total > 0
        n_valid = sum([CASamples.valid]);
    else
        n_valid = 0;
    end
    fprintf('[CA] 收集样本: total=%d  valid=%d\n', n_total, n_valid);
end
