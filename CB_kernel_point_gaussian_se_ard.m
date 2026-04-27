function k = CB_kernel_point_gaussian_se_ard(x, mu, Sigma, kp)
% CB_KERNEL_POINT_GAUSSIAN_SE_ARD  Point-Gaussian ARD SE kernel.

    x = x(:);
    mu = mu(:);
    Sigma = normalize_sigma(Sigma);
    Lambda = diag([kp.ell_x_m^2, kp.ell_y_m^2]);
    A = Lambda + Sigma;
    delta = x - mu;
    k = kp.sigma_f_dB^2 * sqrt(det(Lambda) / max(det(A), realmin)) * ...
        exp(-0.5 * delta' * (A \ delta));
end


function S = normalize_sigma(S)
    if isscalar(S)
        S = eye(2) * S;
    elseif isvector(S)
        S = diag(S(:));
    end
    S = 0.5 * (S + S');
end
