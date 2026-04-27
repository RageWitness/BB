function calibration_events = CB_build_calibration_events_from_source_context(SourceContext, Y_dBm_all, SpatialFP, Config)
% CB_BUILD_CALIBRATION_EVENTS_FROM_SOURCE_CONTEXT  Build CB event structs.

    if nargin < 4
        Config = struct();
    end
    cfg = CB_calib_ap_profile_defaults(Config);
    LC = source_label_constants();

    calibration_events = struct([]);
    if ~isfield(SourceContext, 'label_table') || isempty(SourceContext.label_table)
        return;
    end

    [M, B_obs, T] = size(Y_dBm_all);
    n_added = 0;
    for i = 1:numel(SourceContext.label_table)
        lbl = SourceContext.label_table(i);
        if ~(lbl.label == LC.PERSISTENT_CAL || lbl.label == LC.BROADBAND_CAL || lbl.label == LC.OPPORTUNISTIC)
            continue;
        end

        b = lbl.band_id;
        if b < 1 || b > B_obs || b > SpatialFP.B
            continue;
        end

        [ok_power, P_e, power_reason] = resolve_event_power(lbl, b);
        if ~ok_power
            if cfg.reject_power_none
                fprintf('[CB] skip %s: %s\n', lbl.source_uid, power_reason);
            end
            continue;
        end

        sf = max(1, round(lbl.start_frame));
        ef = min(T, round(lbl.end_frame));
        if ef < sf
            continue;
        end

        rss = squeeze(Y_dBm_all(:, b, sf:ef));
        if isvector(rss)
            rss = rss(:);
        end
        rss_by_ap = finite_mean_dim2(rss);
        if numel(rss_by_ap) ~= M
            continue;
        end

        ev = struct();
        ev.source_uid = lbl.source_uid;
        ev.band_id = b;
        ev.label = lbl.label;
        ev.source_type_name = get_source_type_name(lbl);
        ev.start_frame = sf;
        ev.end_frame = ef;
        ev.time_range = [sf, ef];
        ev.location_prior = lbl.location_prior;
        ev.power_prior = lbl.power_prior;
        ev.P_e_dBm = P_e;
        ev.P0_dBm = get_offline_power(SpatialFP, b);
        ev.rss_dBm_by_ap = rss_by_ap(:)';
        ev.rss_dBm_frames = rss;
        ev.metadata = get_field_default(lbl, 'metadata', struct());
        ev.template_key = get_field_default(lbl, 'template_key', '');

        n_added = n_added + 1;
        if n_added == 1
            calibration_events = ev;
        else
            calibration_events(n_added) = ev; %#ok<AGROW>
        end
    end
end


function [ok, P_e, reason] = resolve_event_power(lbl, band_id)
    ok = false;
    P_e = NaN;
    reason = '';
    if ~isfield(lbl, 'power_prior') || ~isstruct(lbl.power_prior) || ~isfield(lbl.power_prior, 'type')
        reason = 'missing_power_prior';
        return;
    end
    pp = lbl.power_prior;
    switch pp.type
        case 'exact'
            if ~isfield(pp, 'value') || isempty(pp.value)
                reason = 'empty_exact_power';
                return;
            end
            if isscalar(pp.value)
                P_e = pp.value;
            elseif band_id <= numel(pp.value)
                P_e = pp.value(band_id);
            else
                reason = 'exact_power_missing_band';
                return;
            end
        case 'exact_by_band'
            if ~isfield(pp, 'value') || band_id > numel(pp.value)
                reason = 'power_not_exact_for_band';
                return;
            end
            P_e = pp.value(band_id);
        otherwise
            reason = sprintf('unsupported_power_prior_%s', pp.type);
            return;
    end
    ok = isfinite(P_e);
    if ~ok
        reason = 'invalid_power_value';
    end
end


function P0 = get_offline_power(SpatialFP, b)
    if isfield(SpatialFP.band(b), 'reference_power_dBm_used')
        P0 = SpatialFP.band(b).reference_power_dBm_used;
    elseif isfield(SpatialFP, 'reference_power_dBm_used')
        P = SpatialFP.reference_power_dBm_used;
        P0 = P(min(b, numel(P)));
    elseif isfield(SpatialFP, 'ref_power_dBm')
        P = SpatialFP.ref_power_dBm;
        P0 = P(min(b, numel(P)));
    else
        error('[CB] missing offline reference power for band %d', b);
    end
end


function y = finite_mean_dim2(X)
    if isvector(X)
        y = X(:);
        return;
    end
    y = nan(size(X, 1), 1);
    for i = 1:size(X, 1)
        v = X(i, :);
        v = v(isfinite(v));
        if ~isempty(v)
            y(i) = mean(v);
        end
    end
end


function name = get_source_type_name(lbl)
    name = '';
    if isfield(lbl, 'metadata') && isfield(lbl.metadata, 'source_type_name')
        name = lbl.metadata.source_type_name;
    elseif isfield(lbl, 'metadata') && isfield(lbl.metadata, 'label_name')
        name = lbl.metadata.label_name;
    end
end


function v = get_field_default(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end
