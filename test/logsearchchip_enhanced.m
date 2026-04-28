% =========================================================
% 二维空间不确定性投影至一维射线 + 全贝叶斯高斯过程推断 (MCMC-GP)
% =========================================================
clear; clc; close all;

%% 1. 全局参数设置 (物理世界与历史数据)
A_true = 5.5; 
B_true = 2.5;

% 设定射线的角度
theta_ray = pi / 6; 
v_ray = [cos(theta_ray); sin(theta_ray)]; 

% 真实远端点落在射线上，距离为 30m
r_far_true = 30; 
pos_true_2d = r_far_true * v_ray; 

% 生成带空间相关性的历史指纹库数据
N_old = 100;
x_old = linspace(1.5, 40, N_old)';
cov_true = 0.3^2 * exp(-0.5 * (pdist2(x_old, x_old)/5).^2); 
y_old_noisy = A_true + B_true * log(x_old) + chol(cov_true + 0.05*eye(N_old))' * randn(N_old, 1);

% 新版本准观测值 (y值极准)
x1_exact = 1;
y1_exact = A_true + B_true * log(x1_exact); 
yfar_exact = A_true + B_true * log(r_far_true); 
sigma_old = 0.5;
sigma_new = 0.01; 

%% 2. 空间降维预处理模块
% ---------------------------------------------------------
% 【场景 A】 2D高斯 -> 1D高斯
mu_2d = [27; 13]; 
Sigma_2d = [16, 5; 5, 9]; 
W = inv(Sigma_2d);
var_r = 1 / (v_ray' * W * v_ray); 
mu_r = var_r * (v_ray' * W * mu_2d); 
sigma_r = sqrt(var_r);

fprintf('--- 2D高斯投影为1D ---\n');
fprintf('射线方向 1D 均值: %.2fm, 标准差: %.2fm\n\n', mu_r, sigma_r);

% ---------------------------------------------------------
% 【场景 B】 2D矩形 -> 1D线段区间
X_min = 20; X_max = 35;
Y_min = 8;  Y_max = 20;
rx1 = X_min / cos(theta_ray); rx2 = X_max / cos(theta_ray);
ry1 = Y_min / sin(theta_ray); ry2 = Y_max / sin(theta_ray);
r_min = max(min(rx1, rx2), min(ry1, ry2));
r_max = min(max(rx1, rx2), max(ry1, ry2));

fprintf('--- 2D矩形投影为1D ---\n');
fprintf('射线截取 1D 均匀区间: [%.2fm, %.2fm]\n\n', r_min, r_max);

%% 3. 运行全贝叶斯 GP 联合采样 (MCMC-GP)
% ---------------------------------------------------------
num_samples = 3000; % MCMC 采样次数
burn_in = 500;      % 预热期丢弃样本数

fprintf('=== 开始 MCMC 采样计算，请稍候... ===\n');

% (1) 运行 NIGP MCMC (高斯投影)
% 变量: [A, B, sigma_f, l]
theta0_nigp = [4.0, 2.0, 0.5, 5.0]; 
lb_nigp = [3, 0, 0.01, 0.1];
ub_nigp = [8, 4, 2.0,  20.0];
prop_std_nigp = [0.1, 0.1, 0.05, 0.5]; % 提议分布步长

log_target_nigp = @(theta) -gp_nmll_nigp(theta, x_old, y_old_noisy, x1_exact, y1_exact, mu_r, yfar_exact, sigma_old, sigma_new, sigma_r);
[samples_nigp, ~] = mh_sampler(log_target_nigp, theta0_nigp, num_samples, prop_std_nigp, lb_nigp, ub_nigp);

% 提取后验均值
valid_nigp = samples_nigp(burn_in+1:end, :);
A_nigp = mean(valid_nigp(:,1)); B_nigp = mean(valid_nigp(:,2));

% (2) 运行 GP-LVM MCMC (均匀投影隐变量寻优)
% 变量: [A, B, r_lat, sigma_f, l]
theta0_lvm = [4.0, 2.0, (r_min+r_max)/2, 0.5, 5.0]; 
lb_lvm = [3, 0, r_min, 0.01, 0.1];
ub_lvm = [8, 4, r_max, 2.0,  20.0];
prop_std_lvm = [0.1, 0.1, 0.5, 0.05, 0.5];

log_target_lvm = @(theta) -gp_nmll_lvm(theta, x_old, y_old_noisy, x1_exact, y1_exact, yfar_exact, sigma_old, sigma_new);
[samples_lvm, ~] = mh_sampler(log_target_lvm, theta0_lvm, num_samples, prop_std_lvm, lb_lvm, ub_lvm);

valid_lvm = samples_lvm(burn_in+1:end, :);
A_lvm = mean(valid_lvm(:,1)); B_lvm = mean(valid_lvm(:,2)); r_lvm = mean(valid_lvm(:,3));

fprintf('=== 最终 MCMC 后验均值结果 ===\n');
fprintf('真实参数: A = %.4f, B = %.4f\n', A_true, B_true);
fprintf('NIGP (高斯): A = %.4f, B = %.4f\n', A_nigp, B_nigp);
fprintf('GP-LVM(均匀): A = %.4f, B = %.4f, 估计真实距离 = %.4fm\n', A_lvm, B_lvm, r_lvm);

%% 4. 可视化
figure('Position', [100, 100, 1200, 450]);

% 子图 1: 2D 几何空间分布图
subplot(1, 3, 1); hold on; grid on;
patch([X_min X_max X_max X_min], [Y_min Y_min Y_max Y_max], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'm');
plot(mu_2d(1), mu_2d(2), 'b*', 'MarkerSize', 8);
plot([0, 40*cos(theta_ray)], [0, 40*sin(theta_ray)], 'k--', 'LineWidth', 1.5);
plot(pos_true_2d(1), pos_true_2d(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot([r_min*cos(theta_ray), r_max*cos(theta_ray)], [r_min*sin(theta_ray), r_max*sin(theta_ray)], 'm-', 'LineWidth', 4);
plot(mu_r*cos(theta_ray), mu_r*sin(theta_ray), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
axis equal; xlim([0, 40]); ylim([0, 25]);
xlabel('X 坐标 (m)'); ylabel('Y 坐标 (m)'); title('二维几何切片');

% 子图 2: 1D 信号强度拟合图
subplot(1, 3, 2); hold on; grid on;
x_plot = linspace(0.5, 45, 200)';
plot(x_plot, A_true + B_true*log(x_plot), 'k-', 'LineWidth', 2);
scatter(x_old, y_old_noisy, 10, [0.7 0.7 0.7], 'filled');
plot(x_plot, A_nigp + B_nigp*log(x_plot), 'b--', 'LineWidth', 2);
plot(x_plot, A_lvm + B_lvm*log(x_plot), 'm-.', 'LineWidth', 2);
plot(r_lvm, yfar_exact, 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
xlabel('一维距离 r (m)'); ylabel('信号强度 y'); title('MCMC-GP 回归拟合');

% 子图 3: GP-LVM 后验概率分布直方图 (展现参数不确定性)
subplot(1, 3, 3); hold on; grid on;
histogram(valid_lvm(:,2), 30, 'Normalization', 'pdf', 'FaceColor', 'g', 'FaceAlpha', 0.5);
xline(B_true, 'k-', 'LineWidth', 2, 'Label', 'B 真值');
xline(B_lvm, 'r--', 'LineWidth', 2, 'Label', 'B 估计均值');
xlabel('衰减因子 B'); ylabel('后验概率密度'); title('隐变量模型参数 B 后验分布');

%% 辅助函数: Metropolis-Hastings MCMC 采样器
function [samples, accept_rate] = mh_sampler(log_target, theta0, num_samples, prop_std, lb, ub)
    D = length(theta0);
    samples = zeros(num_samples, D);
    theta_curr = theta0;
    logP_curr = log_target(theta_curr);
    accept_count = 0;
    
    for i = 1:num_samples
        % 高斯随机游走提议
        theta_prop = theta_curr + prop_std .* randn(1, D);
        
        % 检查硬边界先验
        if any(theta_prop < lb) || any(theta_prop > ub)
            logP_prop = -Inf;
        else
            logP_prop = log_target(theta_prop);
        end
        
        % 接受概率计算
        if log(rand()) < (logP_prop - logP_curr)
            theta_curr = theta_prop;
            logP_curr = logP_prop;
            accept_count = accept_count + 1;
        end
        samples(i, :) = theta_curr;
    end
    accept_rate = accept_count / num_samples;
end

%% GP 似然计算函数 (计算负似然)
function nmll = gp_nmll_nigp(theta, X_old, Y_old, x1, y1, x_far, y_far, sig_old, sig_new, sig_x)
    A = theta(1); B = theta(2); sig_f = theta(3); l = theta(4);
    X_all = [X_old; x1; x_far]; Y_all = [Y_old; y1; y_far]; N = length(Y_all);
    mu = A + B * log(X_all); res = Y_all - mu;
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    var_far_inflated = sig_new^2 + (B / x_far)^2 * sig_x^2;
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; var_far_inflated]);
    try
        L_chol = chol(K + Noise, 'lower'); alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch, nmll = Inf; end
end

function nmll = gp_nmll_lvm(theta, X_old, Y_old, x1, y1, y_far, sig_old, sig_new)
    A = theta(1); B = theta(2); x_far = theta(3); sig_f = theta(4); l = theta(5);
    X_all = [X_old; x1; x_far]; Y_all = [Y_old; y1; y_far]; N = length(Y_all);
    mu = A + B * log(X_all); res = Y_all - mu;
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; sig_new^2]);
    try
        L_chol = chol(K + Noise, 'lower'); alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch, nmll = Inf; end
end