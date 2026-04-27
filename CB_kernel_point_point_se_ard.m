function k = CB_kernel_point_point_se_ard(x1, x2, kp)
% CB_KERNEL_POINT_POINT_SE_ARD  Point-point ARD SE covariance.

    x1 = x1(:)';
    x2 = x2(:)';
    ell = [kp.ell_x_m, kp.ell_y_m];
    d = (x1 - x2) ./ ell;
    k = kp.sigma_f_dB^2 * exp(-0.5 * sum(d.^2));
end
