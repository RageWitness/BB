function [mu, var_pred] = CA_gpr_predict_se(model, X_pred, ~)
% CA_GPR_PREDICT_SE  用已训练 SE-ARD 模型预测
%
%   model 需包含: X_train, y_train, theta, L, alpha, jitter
%   X_pred: N_pred x 2
%   mu:       N_pred x 1
%   var_pred: N_pred x 1

    ell_x = exp(model.theta(1));
    ell_y = exp(model.theta(2));
    sf    = exp(model.theta(3));
    sn    = exp(model.theta(4));

    X_tr = model.X_train;
    N_pred = size(X_pred, 1);

    % K_star: N_tr x N_pred
    Dx = (X_tr(:,1) - X_pred(:,1)').^2 / ell_x^2;
    Dy = (X_tr(:,2) - X_pred(:,2)').^2 / ell_y^2;
    K_star = sf^2 * exp(-0.5 * (Dx + Dy));

    mu = K_star' * model.alpha;

    % predictive variance
    v = model.L \ K_star;  % L^{-1} K_star
    k_self = sf^2;         % diagonal of K(X_pred, X_pred)
    var_pred = max(0, k_self - sum(v.^2, 1)') + sn^2;
end
