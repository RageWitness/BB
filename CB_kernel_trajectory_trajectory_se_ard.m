function k = CB_kernel_trajectory_trajectory_se_ard(qtraj1, qtraj2, kp)
% CB_KERNEL_TRAJECTORY_TRAJECTORY_SE_ARD  Trajectory-trajectory ARD SE kernel.

    [p1, w1] = CB_trajectory_quadrature(qtraj1, kp.trajectory_quad_n);
    [p2, w2] = CB_trajectory_quadrature(qtraj2, kp.trajectory_quad_n);
    k = 0;
    for i = 1:size(p1, 1)
        for j = 1:size(p2, 1)
            k = k + w1(i) * w2(j) * CB_kernel_point_point_se_ard(p1(i, :), p2(j, :), kp);
        end
    end
end
