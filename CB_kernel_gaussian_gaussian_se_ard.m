function k = CB_kernel_gaussian_gaussian_se_ard(mu1, Sigma1, mu2, Sigma2, kp)
% CB_KERNEL_GAUSSIAN_GAUSSIAN_SE_ARD  Gaussian-Gaussian ARD SE kernel.

    mu1 = mu1(:);
    mu2 = mu2(:);
    Sigma1 = normalize_sigma(Sigma1);
    Sigma2 = normalize_sigma(Sigma2);
    Lambda = diag([kp.ell_x_m^2, kp.ell_y_m^2]);
    A = Lambda + Sigma1 + Sigma2;
    delta = mu1 - mu2;
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
