% =========================================================
% NIGP / GP-LVM: 考虑输入噪声的高斯过程回归参数估计
% 显式均值函数: m(x) = A + B*log(x)
% =========================================================
clear; clc; close all;

%% 1. 全局参数与带空间相关性的数据集生成
A_true = 5.5; 
B_true = 2.5;
xfar_true = 30; 

N_old = 100;
x_old = linspace(1.5, 40, N_old)';
% 生成带局部空间相关性的历史数据 (真实世界信号通常如此)
cov_true = 0.3^2 * exp(-0.5 * (pdist2(x_old, x_old)/5).^2); 
y_old_noisy = A_true + B_true * log(x_old) + chol(cov_true + 0.05*eye(N_old))' * randn(N_old, 1);

% 新版本准观测值
x1_exact = 1;
y1_exact = A_true + B_true * log(x1_exact); 
yfar_exact = A_true + B_true * log(xfar_true); 

%% 2. 场景 A: 高斯分布 (NIGP 误差传播)
fprintf('=== NIGP 高斯分布求解 ===\n');
mu_x = 33;    
sigma_x = 5;  
sigma_old = 0.5;
sigma_new = 0.01; % 准点的基础极小噪声

% 优化变量: theta = [A, B, sigma_f, length_scale] (x_far直接用mu_x，通过方差膨胀处理)
% sigma_f 和 length_scale 是GP核函数的超参数
theta0_nigp = [4.0, 2.0, 0.5, 5.0]; 
lb_nigp = [3, 0, 0.01, 0.1];
ub_nigp = [8, 3, 2.0,  20.0];

% 目标函数: 负边缘对数似然 (NMLL)
obj_nigp = @(theta) gp_nmll_nigp(theta, x_old, y_old_noisy, x1_exact, y1_exact, mu_x, yfar_exact, sigma_old, sigma_new, sigma_x);

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
opt_nigp = fmincon(obj_nigp, theta0_nigp, [], [], [], [], lb_nigp, ub_nigp, [], options);

A_nigp = opt_nigp(1); B_nigp = opt_nigp(2);
fprintf('A = %.4f, B = %.4f (采用方差膨胀吸收位置误差)\n\n', A_nigp, B_nigp);

%% 3. 场景 B: 均匀分布 (GP-LVM 隐变量优化)
fprintf('=== GP-LVM 均匀分布求解 ===\n');
x_a = 21; x_b = 39; 

% 优化变量: theta = [A, B, x_far, sigma_f, length_scale]
theta0_lvm = [4.0, 2.0, 30.0, 0.5, 5.0]; 
lb_lvm = [3, 0, x_a, 0.01, 0.1];
ub_lvm = [8, 3, x_b, 2.0,  20.0];

% 目标函数: NMLL (x_far作为隐变量游走)
obj_lvm = @(theta) gp_nmll_lvm(theta, x_old, y_old_noisy, x1_exact, y1_exact, yfar_exact, sigma_old, sigma_new);

opt_lvm = fmincon(obj_lvm, theta0_lvm, [], [], [], [], lb_lvm, ub_lvm, [], options);

A_lvm = opt_lvm(1); B_lvm = opt_lvm(2); xfar_lvm = opt_lvm(3);
fprintf('A = %.4f, B = %.4f, x_far = %.4f\n\n', A_lvm, B_lvm, xfar_lvm);

%% 4. 绘图对比
x_plot = linspace(0.5, 45, 200)';
y_true = A_true + B_true * log(x_plot);
y_nigp = A_nigp + B_nigp * log(x_plot);
y_lvm  = A_lvm + B_lvm * log(x_plot);

figure('Position', [100, 100, 800, 500]); hold on; grid on;
scatter(x_old, y_old_noisy, 15, [0.7 0.7 0.7], 'filled', 'DisplayName', '带空间相关性的历史数据');
plot(x_plot, y_true, 'k-', 'LineWidth', 2, 'DisplayName', '真实趋势');
plot(x_plot, y_nigp, 'b--', 'LineWidth', 2, 'DisplayName', 'NIGP 均值预测 (高斯)');
plot(x_plot, y_lvm, 'm-.', 'LineWidth', 2, 'DisplayName', 'GP-LVM 均值预测 (均匀)');
plot(xfar_lvm, yfar_exact, 'ms', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'GP-LVM 估计远端位置');
xlabel('距离 x'); ylabel('信号值 y');
title('高斯过程框架下的参数与隐变量估计'); legend('Location','best');

%% 辅助函数: NIGP 的负边缘对数似然 (包含方差膨胀)
function nmll = gp_nmll_nigp(theta, X_old, Y_old, x1, y1, x_far, y_far, sig_old, sig_new, sig_x)
    A = theta(1); B = theta(2); sig_f = theta(3); l = theta(4);
    
    X_all = [X_old; x1; x_far];
    Y_all = [Y_old; y1; y_far];
    N = length(Y_all);
    
    % 计算均值函数残差
    mu = A + B * log(X_all);
    res = Y_all - mu;
    
    % 构建协方差矩阵 (Squared Exponential)
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    
    % NIGP 核心：计算远端点的膨胀方差
    var_far_inflated = sig_new^2 + (B / x_far)^2 * sig_x^2;
    
    % 构建噪声对角阵
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; var_far_inflated]);
    
    % 求解 NMLL (使用Cholesky分解保证数值稳定)
    C = K + Noise;
    try
        L_chol = chol(C, 'lower');
        alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch
        nmll = Inf; % 矩阵非正定则拒绝该参数
    end
end

%% 辅助函数: GP-LVM 的负边缘对数似然 (x_far 作为未知参数)
function nmll = gp_nmll_lvm(theta, X_old, Y_old, x1, y1, y_far, sig_old, sig_new)
    A = theta(1); B = theta(2); x_far = theta(3); sig_f = theta(4); l = theta(5);
    
    X_all = [X_old; x1; x_far];
    Y_all = [Y_old; y1; y_far];
    N = length(Y_all);
    
    mu = A + B * log(X_all);
    res = Y_all - mu;
    
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; sig_new^2]);
    
    C = K + Noise;
    try
        L_chol = chol(C, 'lower');
        alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch
        nmll = Inf;
    end
end