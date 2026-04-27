function k = CB_kernel_rect_trajectory_se_ard(bbox, qtraj, kp)
% CB_KERNEL_RECT_TRAJECTORY_SE_ARD  Rectangle-trajectory ARD SE kernel.

    [pts, w] = CB_trajectory_quadrature(qtraj, kp.trajectory_quad_n);
    k = 0;
    for i = 1:size(pts, 1)
        k = k + w(i) * CB_kernel_point_rect_se_ard(pts(i, :), bbox, kp);
    end
end
