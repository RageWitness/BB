function k = CB_kernel_rect_gaussian_se_ard(bbox, mu, Sigma, kp)
% CB_KERNEL_RECT_GAUSSIAN_SE_ARD  Rectangle-Gaussian ARD SE kernel.

    bbox = bbox(:)';
    mu = mu(:)';
    Sigma = normalize_sigma(Sigma);

    offdiag = Sigma - diag(diag(Sigma));
    if max(abs(offdiag(:))) < 1e-12
        a = [bbox(1), bbox(3)];
        b = [bbox(2), bbox(4)];
        L = b - a;
        ell = [kp.ell_x_m, kp.ell_y_m];
        sig2 = diag(Sigma)';
        if any(L <= max(kp.eps_dist, 1e-9))
            xc = 0.5 * (a + b);
            k = CB_kernel_point_gaussian_se_ard(xc, mu, Sigma, kp);
            return;
        end

        prod_term = 1;
        for r = 1:2
            den = sqrt(2 * (ell(r)^2 + sig2(r)));
            term = sqrt(pi / 2) * ell(r) / L(r) * ...
                (erf((b(r) - mu(r)) / den) - erf((a(r) - mu(r)) / den));
            prod_term = prod_term * term;
        end
        k = kp.sigma_f_dB^2 * prod_term;
    else
        k = rect_gaussian_quadrature(bbox, mu, Sigma, kp);
    end
end


function k = rect_gaussian_quadrature(bbox, mu, Sigma, kp)
    n = max(3, kp.rect_quad_n);
    xs = linspace(bbox(1), bbox(2), n);
    ys = linspace(bbox(3), bbox(4), n);
    wx = trapezoid_weights(n);
    wy = trapezoid_weights(n);
    val = 0;
    wsum = 0;
    for ix = 1:n
        for iy = 1:n
            x = [xs(ix), ys(iy)];
            w = wx(ix) * wy(iy);
            val = val + w * CB_kernel_point_gaussian_se_ard(x, mu, Sigma, kp);
            wsum = wsum + w;
        end
    end
    k = val / max(wsum, realmin);
end


function w = trapezoid_weights(n)
    w = ones(1, n);
    w([1, end]) = 0.5;
end


function S = normalize_sigma(S)
    if isscalar(S)
        S = eye(2) * S;
    elseif isvector(S)
        S = diag(S(:));
    end
    S = 0.5 * (S + S');
end
