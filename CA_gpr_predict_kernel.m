function [mu, var_pred] = CA_gpr_predict_kernel(model, X_pred, KernelCtx, Config)
% CA_GPR_PREDICT_KERNEL  Predict with the configured CA kernel.

    X_pred = ensure_2d(X_pred);
    if isempty(X_pred)
        mu = [];
        var_pred = [];
        return;
    end

    K_star = CA_kernel_matrix( ...
        X_pred, model.X_train, model.theta, ...
        model.band_id, model.ap_id, KernelCtx, Config);  % N_pred x N_train

    mu = K_star * model.alpha;

    v = model.L \ K_star';
    sf = exp(model.theta(3));
    var_pred = sf^2 - sum(v.^2, 1)';
    var_pred = max(var_pred, 0);
end


function X = ensure_2d(X)
    if isempty(X)
        X = zeros(0, 2);
    elseif size(X, 2) ~= 2 && size(X, 1) == 2
        X = X';
    end
end
