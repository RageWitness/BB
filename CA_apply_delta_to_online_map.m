function [OnlineSpatialFP_new, UpdateLog] = CA_apply_delta_to_online_map(OnlineSpatialFP, DeltaA, Config)
% CA_APPLY_DELTA_TO_ONLINE_MAP  回写 centered_dBm 并记录 UpdateLog
%
%   Z_new(b,m,g) = Z_old(b,m,g) + lambda(b,m,g) * clip(mu(b,m,g))

    Config = CA_fill_defaults(Config);
    max_delta = Config.ca.update.max_delta_dB;

    OnlineSpatialFP_new = OnlineSpatialFP;
    B = numel(DeltaA.band);

    total_updated = 0;
    all_applied = [];

    for b = 1:B
        mu_raw  = DeltaA.band(b).mu;
        lam     = DeltaA.band(b).lambda;
        valid   = DeltaA.band(b).valid;

        clipped_mu = max(-max_delta, min(max_delta, mu_raw));
        applied = lam .* clipped_mu;
        applied(~valid) = 0;

        OnlineSpatialFP_new.band(b).centered_dBm = ...
            OnlineSpatialFP.band(b).centered_dBm + applied;

        n_b = sum(valid(:));
        total_updated = total_updated + n_b;

        vals = applied(valid);
        if ~isempty(vals)
            all_applied = [all_applied; vals(:)]; %#ok<AGROW>
        end
    end

    UpdateLog = struct();
    UpdateLog.n_updated_entries = total_updated;

    if total_updated > 0
        UpdateLog.mean_abs_delta = mean(abs(all_applied));
        UpdateLog.max_abs_delta  = max(abs(all_applied));
    else
        UpdateLog.mean_abs_delta = 0;
        UpdateLog.max_abs_delta  = 0;
    end
    UpdateLog.status = 'ok';
    UpdateLog.timestamp = now;

    OnlineSpatialFP_new.update_round = OnlineSpatialFP.update_round + 1;
    if isempty(OnlineSpatialFP_new.update_history)
        OnlineSpatialFP_new.update_history = UpdateLog;
    else
        OnlineSpatialFP_new.update_history(end+1) = UpdateLog;
    end

    fprintf('[CA] 回写完成: updated=%d  mean_abs=%.4f  max_abs=%.4f\n', ...
        total_updated, UpdateLog.mean_abs_delta, UpdateLog.max_abs_delta);
end
