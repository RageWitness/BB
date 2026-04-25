function result = CA_calibrate_fp_by_trusted_source( ...
    ap_xy, grid_xy, source_xy, rss_dB_frames, F_cur, power_mode, cfg)
% CA_CALIBRATE_FP_BY_TRUSTED_SOURCE  Single trusted-source RF_raw calibration.

    if nargin < 7 || isempty(cfg)
        cfg = CA_trusted_rss_defaults();
    else
        cfg = fill_missing_cfg(cfg);
    end

    result = init_result(power_mode, F_cur);

    search_result = CA_estimate_n0_top2_ratio_search( ...
        ap_xy, source_xy, rss_dB_frames, cfg);
    result.search_result = search_result;

    if ~search_result.valid
        result.warning = search_result.warning;
        return;
    end

    F_tmp = CA_build_tmp_fp_from_trusted_source( ...
        ap_xy, grid_xy, source_xy, search_result.Rbar_dB, ...
        search_result.n_hat, cfg);

    result.F_tmp = F_tmp;
    result.n_hat = search_result.n_hat;

    switch lower(power_mode)
        case 'known_same'
            delta_pow = 0;
            gamma_hat = NaN;

        case 'known_different'
            require_field(cfg, 'P_c');
            require_field(cfg, 'P_off');
            delta_pow = cfg.P_off - cfg.P_c;
            gamma_hat = NaN;

        case 'unknown'
            delta_pow = median_omitnan(F_cur(:) - F_tmp(:));
            gamma_hat = delta_pow;
            if ~isfinite(delta_pow)
                result.warning = 'invalid_unknown_power_shift';
                return;
            end

        otherwise
            result.warning = sprintf('unknown_power_mode_%s', power_mode);
            return;
    end

    result.delta_pow = delta_pow;
    result.gamma_hat = gamma_hat;
    result.F_new = F_tmp + delta_pow;
    result.valid = true;
    result.warning = search_result.warning;
end


function cfg = fill_missing_cfg(cfg)
    d = CA_trusted_rss_defaults();
    names = fieldnames(d);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k}) || isempty(cfg.(names{k}))
            cfg.(names{k}) = d.(names{k});
        end
    end
end


function result = init_result(power_mode, F_cur)
    result = struct();
    result.F_new = F_cur;
    result.F_tmp = [];
    result.n_hat = NaN;
    result.delta_pow = NaN;
    result.gamma_hat = NaN;
    result.power_mode = power_mode;
    result.valid = false;
    result.warning = '';
    result.search_result = struct();
end


function require_field(cfg, fname)
    if ~isfield(cfg, fname) || isempty(cfg.(fname)) || ...
            ~isscalar(cfg.(fname)) || ~isfinite(cfg.(fname))
        error('[CA] cfg.%s must be a finite scalar', fname);
    end
end


function y = median_omitnan(x)
    x = x(isfinite(x));
    if isempty(x)
        y = NaN;
    else
        y = median(x);
    end
end
