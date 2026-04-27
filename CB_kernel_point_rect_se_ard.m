function k = CB_kernel_point_rect_se_ard(x, bbox, kp)
% CB_KERNEL_POINT_RECT_SE_ARD  Point-rectangle uncertain-input kernel.

    x = x(:)';
    bbox = bbox(:)';
    a = [bbox(1), bbox(3)];
    b = [bbox(2), bbox(4)];
    ell = [kp.ell_x_m, kp.ell_y_m];

    if any((b - a) <= max(kp.eps_dist, 1e-9))
        xc = 0.5 * (a + b);
        k = CB_kernel_point_point_se_ard(x, xc, kp);
        return;
    end

    prod_term = 1;
    for r = 1:2
        Ar_b = A_fun(b(r) - x(r), ell(r));
        Ar_a = A_fun(a(r) - x(r), ell(r));
        prod_term = prod_term * ((Ar_b - Ar_a) / (b(r) - a(r)));
    end
    k = kp.sigma_f_dB^2 * prod_term;
end


function y = A_fun(z, ell)
    y = sqrt(pi / 2) * ell * erf(z / (sqrt(2) * ell));
end
