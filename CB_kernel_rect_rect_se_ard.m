function k = CB_kernel_rect_rect_se_ard(bbox1, bbox2, kp)
% CB_KERNEL_RECT_RECT_SE_ARD  Rectangle-rectangle ARD SE integral kernel.

    bbox1 = bbox1(:)';
    bbox2 = bbox2(:)';
    ai = [bbox1(1), bbox1(3)];
    bi = [bbox1(2), bbox1(4)];
    aj = [bbox2(1), bbox2(3)];
    bj = [bbox2(2), bbox2(4)];
    Li = bi - ai;
    Lj = bj - aj;
    ell = [kp.ell_x_m, kp.ell_y_m];

    if any(Li <= max(kp.eps_dist, 1e-9)) || any(Lj <= max(kp.eps_dist, 1e-9))
        xi = 0.5 * (ai + bi);
        xj = 0.5 * (aj + bj);
        k = CB_kernel_point_point_se_ard(xi, xj, kp);
        return;
    end

    prod_term = 1;
    for r = 1:2
        num = B_fun(bi(r) - aj(r), ell(r)) ...
            - B_fun(ai(r) - aj(r), ell(r)) ...
            - B_fun(bi(r) - bj(r), ell(r)) ...
            + B_fun(ai(r) - bj(r), ell(r));
        prod_term = prod_term * (num / (Li(r) * Lj(r)));
    end
    k = kp.sigma_f_dB^2 * prod_term;
end


function y = B_fun(z, ell)
    y = z * sqrt(pi / 2) * ell * erf(z / (sqrt(2) * ell)) ...
        + ell^2 * exp(-(z.^2) / (2 * ell^2));
end
