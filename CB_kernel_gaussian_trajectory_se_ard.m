function k = CB_kernel_gaussian_trajectory_se_ard(mu, Sigma, qtraj, kp)
% CB_KERNEL_GAUSSIAN_TRAJECTORY_SE_ARD  Gaussian-trajectory ARD SE kernel.

    [pts, w] = CB_trajectory_quadrature(qtraj, kp.trajectory_quad_n);
    k = 0;
    for i = 1:size(pts, 1)
        k = k + w(i) * CB_kernel_point_gaussian_se_ard(pts(i, :), mu, Sigma, kp);
    end
end
