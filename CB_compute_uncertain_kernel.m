function kbar = CB_compute_uncertain_kernel(q_i, q_j, kernel_params)
% CB_COMPUTE_UNCERTAIN_KERNEL  Unified uncertain-input ARD SE kernel.

    kp = normalize_kernel_params(kernel_params);
    ti = lower(q_i.type);
    tj = lower(q_j.type);

    switch ti
        case 'point'
            switch tj
                case 'point'
                    kbar = CB_kernel_point_point_se_ard(q_i.x, q_j.x, kp);
                case 'rect'
                    kbar = CB_kernel_point_rect_se_ard(q_i.x, q_j.bbox, kp);
                case 'gaussian'
                    kbar = CB_kernel_point_gaussian_se_ard(q_i.x, q_j.mu, q_j.Sigma, kp);
                case 'trajectory'
                    kbar = CB_kernel_point_trajectory_se_ard(q_i.x, q_j, kp);
                otherwise
                    error('[CB] unsupported q type: %s', tj);
            end

        case 'rect'
            switch tj
                case 'point'
                    kbar = CB_kernel_point_rect_se_ard(q_j.x, q_i.bbox, kp);
                case 'rect'
                    kbar = CB_kernel_rect_rect_se_ard(q_i.bbox, q_j.bbox, kp);
                case 'gaussian'
                    kbar = CB_kernel_rect_gaussian_se_ard(q_i.bbox, q_j.mu, q_j.Sigma, kp);
                case 'trajectory'
                    kbar = CB_kernel_rect_trajectory_se_ard(q_i.bbox, q_j, kp);
                otherwise
                    error('[CB] unsupported q type: %s', tj);
            end

        case 'gaussian'
            switch tj
                case 'point'
                    kbar = CB_kernel_point_gaussian_se_ard(q_j.x, q_i.mu, q_i.Sigma, kp);
                case 'rect'
                    kbar = CB_kernel_rect_gaussian_se_ard(q_j.bbox, q_i.mu, q_i.Sigma, kp);
                case 'gaussian'
                    kbar = CB_kernel_gaussian_gaussian_se_ard(q_i.mu, q_i.Sigma, q_j.mu, q_j.Sigma, kp);
                case 'trajectory'
                    kbar = CB_kernel_gaussian_trajectory_se_ard(q_i.mu, q_i.Sigma, q_j, kp);
                otherwise
                    error('[CB] unsupported q type: %s', tj);
            end

        case 'trajectory'
            switch tj
                case 'point'
                    kbar = CB_kernel_point_trajectory_se_ard(q_j.x, q_i, kp);
                case 'rect'
                    kbar = CB_kernel_rect_trajectory_se_ard(q_j.bbox, q_i, kp);
                case 'gaussian'
                    kbar = CB_kernel_gaussian_trajectory_se_ard(q_j.mu, q_j.Sigma, q_i, kp);
                case 'trajectory'
                    kbar = CB_kernel_trajectory_trajectory_se_ard(q_i, q_j, kp);
                otherwise
                    error('[CB] unsupported q type: %s', tj);
            end

        otherwise
            error('[CB] unsupported q type: %s', ti);
    end

    if ~isfinite(kbar)
        kbar = 0;
    end
end


function kp = normalize_kernel_params(kp)
    kp = set_default(kp, 'sigma_f_dB', 8);
    kp = set_default(kp, 'ell_x_m', 35);
    kp = set_default(kp, 'ell_y_m', 35);
    kp = set_default(kp, 'eps_dist', 1e-6);
    kp = set_default(kp, 'trajectory_quad_n', 25);
    kp = set_default(kp, 'rect_quad_n', 9);
end


function s = set_default(s, name, value)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = value;
    end
end
