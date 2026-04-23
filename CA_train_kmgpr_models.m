function Models = CA_train_kmgpr_models(CADataset, Partition, Config)
% CA_TRAIN_KMGPR_MODELS  对每个 (b,m,k) 训练 GPR
%
%   Models.band(b).ap(m).cluster(k):
%       .theta, .L, .alpha, .X_train, .y_train
%       .n_train, .status, .valid, .jitter
%       .learned_hyperparams (struct with ell_x, ell_y, sigma_f, sigma_n)

    Config = CA_fill_defaults(Config);
    gpr_cfg = Config.ca.gpr;
    B = numel(CADataset.band);
    K = Partition.K;

    for b = 1:B
        X_all = CADataset.band(b).X;
        Y_all = CADataset.band(b).Y;   % N x M
        cid   = CADataset.band(b).cluster_id;
        N = size(X_all, 1);

        if N == 0
            M_ap = 0;
        else
            M_ap = size(Y_all, 2);
        end

        for m = 1:M_ap
            for k = 1:K
                mdl = struct();
                mdl.valid = false;
                mdl.status = '';
                mdl.n_train = 0;
                mdl.theta = [];
                mdl.L = [];
                mdl.alpha = [];
                mdl.X_train = [];
                mdl.y_train = [];
                mdl.jitter = gpr_cfg.jitter_init;
                mdl.learned_hyperparams = struct();

                % select cluster samples with finite residual for this AP
                if N > 0
                    in_cluster = (cid == k);
                    y_col = Y_all(:, m);
                    valid_rows = in_cluster & isfinite(y_col);
                    X_k = X_all(valid_rows, :);
                    y_k = y_col(valid_rows);
                else
                    X_k = [];
                    y_k = [];
                end

                n_k = numel(y_k);
                mdl.n_train = n_k;

                if n_k < gpr_cfg.n_min
                    mdl.status = 'insufficient_cluster_samples';
                    Models.band(b).ap(m).cluster(k) = mdl;
                    continue;
                end

                % hyperparameter learning
                if gpr_cfg.learn_hyperparams && n_k >= gpr_cfg.n_hyper_min
                    [theta, opt_info] = CA_optimize_gpr_hyperparams(X_k, y_k, Config);
                    mdl.status = opt_info.status;
                else
                    theta = [ log(gpr_cfg.init.length_scale_x), ...
                              log(gpr_cfg.init.length_scale_y), ...
                              log(gpr_cfg.init.sigma_f), ...
                              log(gpr_cfg.init.sigma_n) ];
                    if n_k < gpr_cfg.n_hyper_min
                        mdl.status = 'default_hyperparams_due_to_few_samples';
                    else
                        mdl.status = 'hyperparams_learning_disabled';
                    end
                end

                ell_x = exp(theta(1));
                ell_y = exp(theta(2));
                sf    = exp(theta(3));
                sn    = exp(theta(4));

                mdl.learned_hyperparams.ell_x   = ell_x;
                mdl.learned_hyperparams.ell_y   = ell_y;
                mdl.learned_hyperparams.sigma_f = sf;
                mdl.learned_hyperparams.sigma_n = sn;

                % build kernel and cache L, alpha
                Dx = (X_k(:,1) - X_k(:,1)').^2 / ell_x^2;
                Dy = (X_k(:,2) - X_k(:,2)').^2 / ell_y^2;
                K_mat = sf^2 * exp(-0.5 * (Dx + Dy)) + sn^2 * eye(n_k);

                j = gpr_cfg.jitter_init;
                chol_ok = false;
                while j <= gpr_cfg.jitter_max
                    try
                        L = chol(K_mat + j * eye(n_k), 'lower');
                        chol_ok = true;
                        break;
                    catch
                        j = j * gpr_cfg.jitter_factor;
                    end
                end

                if ~chol_ok
                    mdl.status = 'cholesky_failed';
                    Models.band(b).ap(m).cluster(k) = mdl;
                    continue;
                end

                mdl.theta   = theta;
                mdl.L       = L;
                mdl.alpha   = L' \ (L \ y_k);
                mdl.X_train = X_k;
                mdl.y_train = y_k;
                mdl.jitter  = j;
                mdl.valid   = true;

                Models.band(b).ap(m).cluster(k) = mdl;
            end
        end

        % fill missing M_ap for empty bands
        if M_ap == 0
            Models.band(b).ap = struct('cluster', {});
        end
    end

    % summary
    n_valid = 0; n_total = 0;
    for b = 1:B
        for m = 1:numel(Models.band(b).ap)
            for k = 1:K
                n_total = n_total + 1;
                if Models.band(b).ap(m).cluster(k).valid
                    n_valid = n_valid + 1;
                end
            end
        end
    end
    fprintf('[CA] GPR 模型: total=%d  valid=%d\n', n_total, n_valid);
end
