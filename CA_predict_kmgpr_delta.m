function DeltaA = CA_predict_kmgpr_delta(Models, CADataset, OnlineSpatialFP, Partition, Config)
% CA_PREDICT_KMGPR_DELTA  cluster 模型预测 + 边界平滑融合
%
%   DeltaA.band(b).mu:    M x G
%   DeltaA.band(b).var:   M x G
%   DeltaA.band(b).valid: M x G logical
%   DeltaA.band(b).n_eff: M x G

    Config = CA_fill_defaults(Config);
    B = OnlineSpatialFP.B;
    M_ap = OnlineSpatialFP.M;
    G = OnlineSpatialFP.G;
    grid_xy = OnlineSpatialFP.grid_xy;
    K = Partition.K;
    centers = Partition.centers;

    enable_blend = Config.ca.kmeans.enable_blend;
    blend_k = Config.ca.kmeans.blend_k;
    r_blend = Config.ca.kmeans.blend_radius;
    skip_var = Config.ca.sample.skip_var;
    eps_v = 1e-12;

    % precompute grid-to-center distances
    g2c_dist2 = zeros(G, K);
    for k = 1:K
        diffs = grid_xy - centers(k, :);
        g2c_dist2(:, k) = sum(diffs.^2, 2);
    end

    for b = 1:B
        mu_out    = zeros(M_ap, G);
        var_out   = ones(M_ap, G) * skip_var;
        valid_out = false(M_ap, G);
        n_eff_out = zeros(M_ap, G);

        has_models = b <= numel(Models.band) && numel(Models.band(b).ap) >= M_ap;
        if ~has_models
            DeltaA.band(b).mu     = mu_out;
            DeltaA.band(b).var    = var_out;
            DeltaA.band(b).valid  = valid_out;
            DeltaA.band(b).n_eff  = n_eff_out;
            DeltaA.band(b).status = 'no_valid_samples_for_band';
            DeltaA.band(b).status_per_ap = repmat({'no_valid_samples_for_ap'}, 1, M_ap);
            continue;
        end

        status_per_ap = cell(1, M_ap);

        for m = 1:M_ap
            any_valid = false;
            for k = 1:K
                if Models.band(b).ap(m).cluster(k).valid
                    any_valid = true; break;
                end
            end
            if ~any_valid
                status_per_ap{m} = 'no_valid_samples_for_ap';
                continue;
            end

            status_per_ap{m} = 'ok';

            if enable_blend && blend_k > 1
                % --- blended prediction per grid point ---
                for g = 1:G
                    [~, sorted_k] = sort(g2c_dist2(g, :));
                    cands = sorted_k(1:min(blend_k, K));

                    w_total = 0;
                    mu_acc = 0;
                    inv_var_acc = 0;
                    n_eff_acc = 0;

                    for ci = 1:numel(cands)
                        kk = cands(ci);
                        mdl = Models.band(b).ap(m).cluster(kk);
                        if ~mdl.valid, continue; end

                        [mu_k, var_k] = CA_gpr_predict_se(mdl, grid_xy(g, :), Config);

                        w_s = exp(-g2c_dist2(g, kk) / (2 * r_blend^2));
                        w_v = 1 / (var_k + eps_v);
                        w = w_s * w_v;

                        mu_acc = mu_acc + w * mu_k;
                        inv_var_acc = inv_var_acc + w_v;
                        w_total = w_total + w;
                        if isfield(mdl, 'n_eff') && mdl.n_eff > 0
                            n_eff_acc = n_eff_acc + w_s * mdl.n_eff;
                        else
                            n_eff_acc = n_eff_acc + w_s * mdl.n_train;
                        end
                    end

                    if w_total > 0
                        mu_out(m, g) = mu_acc / w_total;
                        var_out(m, g) = 1 / (inv_var_acc + eps_v);
                        valid_out(m, g) = true;
                        n_eff_out(m, g) = n_eff_acc;
                    end
                end
            else
                % --- hard assignment (no blend) ---
                grid_cid = Partition.cluster_id;
                for k = 1:K
                    mdl = Models.band(b).ap(m).cluster(k);
                    if ~mdl.valid, continue; end
                    g_mask = (grid_cid == k);
                    g_idx = find(g_mask);
                    if isempty(g_idx), continue; end

                    [mu_k, var_k] = CA_gpr_predict_se(mdl, grid_xy(g_idx, :), Config);

                    mu_out(m, g_idx) = mu_k';
                    var_out(m, g_idx) = var_k';
                    valid_out(m, g_idx) = true;
                    if isfield(mdl, 'n_eff') && mdl.n_eff > 0
                        n_eff_out(m, g_idx) = mdl.n_eff;
                    else
                        n_eff_out(m, g_idx) = mdl.n_train;
                    end
                end
            end
        end

        n_valid_aps = sum(~strcmp(status_per_ap, 'no_valid_samples_for_ap'));
        if n_valid_aps == M_ap
            band_status = 'ok';
        elseif n_valid_aps > 0
            band_status = 'partial_ap_update';
        else
            band_status = 'no_valid_samples_for_band';
        end

        DeltaA.band(b).mu     = mu_out;
        DeltaA.band(b).var    = var_out;
        DeltaA.band(b).valid  = valid_out;
        DeltaA.band(b).n_eff  = n_eff_out;
        DeltaA.band(b).status = band_status;
        DeltaA.band(b).status_per_ap = status_per_ap;
    end
end
