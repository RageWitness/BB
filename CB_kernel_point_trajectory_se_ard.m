function k = CB_kernel_point_trajectory_se_ard(x, qtraj, kp)
% CB_KERNEL_POINT_TRAJECTORY_SE_ARD  Point-trajectory ARD SE kernel.

    [pts, w] = CB_trajectory_quadrature(qtraj, kp.trajectory_quad_n);
    k = 0;
    for i = 1:size(pts, 1)
        k = k + w(i) * CB_kernel_point_point_se_ard(x, pts(i, :), kp);
    end
end
