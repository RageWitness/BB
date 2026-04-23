function [z_centered, y_dBm_mean, y_lin_mean, valid_ap_mask, info] = ...
    CA_compute_centered_observation(obs_segment_lin, Config)
% CA_COMPUTE_CENTERED_OBSERVATION  清洗、聚合、centered（per-AP 独立）
%
%   obs_segment_lin: M x T 线性功率矩阵
%   返回 M x 1 向量（无效 AP 为 NaN）

    Config = CA_fill_defaults(Config);
    min_aps = Config.ca.sample.min_valid_aps;

    [M, ~] = size(obs_segment_lin);

    info = struct('status', 'ok', 'n_valid_aps', 0, 'obs_mean_power', NaN);

    invalid_mask = isnan(obs_segment_lin) | isinf(obs_segment_lin) | (obs_segment_lin <= 0);
    obs_clean = obs_segment_lin;
    obs_clean(invalid_mask) = NaN;

    y_lin_mean = zeros(M, 1);
    y_dBm_mean = nan(M, 1);
    for m = 1:M
        row = obs_clean(m, :);
        valid_t = ~isnan(row);
        if any(valid_t)
            y_lin_mean(m) = mean(row(valid_t));
            y_dBm_mean(m) = 10 * log10(y_lin_mean(m));
        else
            y_lin_mean(m) = NaN;
            y_dBm_mean(m) = NaN;
        end
    end

    valid_ap_mask = isfinite(y_dBm_mean);
    n_valid = sum(valid_ap_mask);
    info.n_valid_aps = n_valid;

    z_centered = nan(M, 1);

    if n_valid < min_aps
        info.status = 'too_few_valid_aps';
        return;
    end

    mu_r = mean(y_dBm_mean(valid_ap_mask));
    info.obs_mean_power = mu_r;

    z_centered(valid_ap_mask) = y_dBm_mean(valid_ap_mask) - mu_r;
end
