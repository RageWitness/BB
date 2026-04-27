function mu = CB_eval_ap_profile_mean(mean_model, x)
% CB_EVAL_AP_PROFILE_MEAN  Evaluate log-distance AP-profile mean.

    if isempty(x)
        mu = [];
        return;
    end
    if size(x, 2) ~= 2
        x = x';
    end

    if isfield(mean_model, 'ap_xy') && all(isfinite(mean_model.ap_xy))
        d = sqrt(sum((x - mean_model.ap_xy).^2, 2));
        d = max(d, 1.0);
        mu = mean_model.intercept_dB + mean_model.slope_dB_per_log10m * log10(d);
    else
        mu = mean_model.intercept_dB * ones(size(x, 1), 1);
    end
end
