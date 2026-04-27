function report = CB_validate_ap_profile_calibration()
% CB_VALIDATE_AP_PROFILE_CALIBRATION  Basic CB kernel and update checks.

    report = struct();
    report.tests = {};
    report.pass = true;

    cfg = CB_calib_ap_profile_defaults(struct());
    kp = struct('sigma_f_dB', 5, 'ell_x_m', 20, 'ell_y_m', 30, ...
        'eps_dist', 1e-9, 'trajectory_quad_n', 15, 'rect_quad_n', 7);

    try
        q1 = struct('type', 'point', 'x', [10, 20]);
        q2 = struct('type', 'point', 'x', [15, 23]);
        qrect0 = struct('type', 'rect', 'bbox', [10, 10+1e-9, 20, 20+1e-9]);
        qg0 = struct('type', 'gaussian', 'mu', [10, 20], 'Sigma', 1e-12 * eye(2));

        kpp = CB_compute_uncertain_kernel(q1, q1, kp);
        kpr0 = CB_compute_uncertain_kernel(q1, qrect0, kp);
        krr0 = CB_compute_uncertain_kernel(qrect0, qrect0, kp);
        kpg0 = CB_compute_uncertain_kernel(q1, qg0, kp);
        kgg0 = CB_compute_uncertain_kernel(qg0, qg0, kp);

        add_test(abs(kpr0 - kpp) < 1e-6, 'point-rect zero-area degenerates');
        add_test(abs(krr0 - kpp) < 1e-6, 'rect-rect zero-area degenerates');
        add_test(abs(kpg0 - kpp) < 1e-6, 'point-gaussian zero-cov degenerates');
        add_test(abs(kgg0 - kpp) < 1e-6, 'gaussian-gaussian zero-cov degenerates');

        kab = CB_compute_uncertain_kernel(q1, q2, kp);
        kba = CB_compute_uncertain_kernel(q2, q1, kp);
        add_test(abs(kab - kba) < 1e-10, 'kernel symmetry point-point');

        qs = {q1, q2, qrect0, qg0};
        K = zeros(numel(qs));
        for i = 1:numel(qs)
            for j = 1:numel(qs)
                K(i,j) = CB_compute_uncertain_kernel(qs{i}, qs{j}, kp);
            end
        end
        K = 0.5 * (K + K') + 1e-8 * eye(numel(qs));
        [~, p] = chol(K);
        add_test(p == 0, 'kernel covariance PSD with jitter');

        y_same = normalize_power(-55, 20, 20);
        y_diff = normalize_power(-55, 30, 20);
        add_test(abs(y_same + 55) < 1e-12, 'power normalization P_e=P0');
        add_test(abs(y_diff + 65) < 1e-12, 'power normalization P_e higher by 10 dB');

        add_test(cfg.power_known_only, 'power-known-only default enabled');
    catch ME
        report.pass = false;
        report.error = ME.message;
    end

    fprintf('[CB-Validate] pass=%d, tests=%d\n', report.pass, numel(report.tests));

    function add_test(ok, name)
        item = struct('name', name, 'pass', logical(ok));
        report.tests{end+1} = item;
        if ~ok
            report.pass = false;
            fprintf('[CB-Validate] FAIL: %s\n', name);
        else
            fprintf('[CB-Validate] PASS: %s\n', name);
        end
    end
end


function ytilde = normalize_power(y, Pe, P0)
    ytilde = y - (Pe - P0);
end
