% =========================================================
% 二维空间不确定性投影至一维射线 + NIGP/GP-LVM 曲线拟合
% =========================================================
clear; clc; close all;

%% 1. 全局参数设置 (物理世界与历史数据)
A_true = 5.5; 
B_true = 2.1;

% 设定射线的角度 (假设参考基准位于原点，我们沿着30度方向测量)
theta_ray = pi / 6; 
v_ray = [cos(theta_ray); sin(theta_ray)]; % 射线方向向量

% 真实远端点落在射线上，假设距离为 30m
r_far_true = 30; 
pos_true_2d = r_far_true * v_ray; % 真实的二维坐标

% 生成历史指纹库数据 (带空间相关性)
N_old = 100;
x_old = linspace(1.5, 40, N_old)';
cov_true = 0.3^2 * exp(-0.5 * (pdist2(x_old, x_old)/5).^2); 
y_old_noisy = A_true + B_true * log(x_old) + chol(cov_true + 0.05*eye(N_old))' * randn(N_old, 1);

% 新版本准观测值 (y值极准)
x1_exact = 1;
y1_exact = A_true + B_true * log(x1_exact); 
yfar_exact = A_true + B_true * log(r_far_true); 
sigma_old = 0.5;
sigma_new = 0.01; % 新观测值的极小噪声

%% 2. 空间降维预处理模块
% ---------------------------------------------------------
% 【场景 A】 二维高斯分布 -> 射线投影 -> 一维高斯
% 假设测量系统的二维误差为:
mu_2d = [27; 13]; % 二维高斯中心
Sigma_2d = [16, 5; 5, 9]; % 二维协方差矩阵

% 数学投影：计算射线上的一维条件高斯分布参数
W = inv(Sigma_2d);
var_r = 1 / (v_ray' * W * v_ray); % 一维射线上方差
mu_r = var_r * (v_ray' * W * mu_2d); % 一维射线上的期望中心
sigma_r = sqrt(var_r);

fprintf('--- 2D高斯投影为1D ---\n');
fprintf('射线方向 1D 均值: %.2fm, 标准差: %.2fm\n\n', mu_r, sigma_r);

% ---------------------------------------------------------
% 【场景 B】 二维矩形均匀区域 -> 射线相交 -> 一维线段区间
% 假设目标分布在一个矩形内: X在[20, 35], Y在[8, 20]
X_min = 20; X_max = 35;
Y_min = 8;  Y_max = 20;

% 计算射线与X边界交点的距离
rx1 = X_min / cos(theta_ray); rx2 = X_max / cos(theta_ray);
% 计算射线与Y边界交点的距离
ry1 = Y_min / sin(theta_ray); ry2 = Y_max / sin(theta_ray);

% 矩形与射线相交的有效线段(求交集)
r_min = max(min(rx1, rx2), min(ry1, ry2));
r_max = min(max(rx1, rx2), max(ry1, ry2));

fprintf('--- 2D矩形投影为1D ---\n');
fprintf('射线截取 1D 均匀区间: [%.2fm, %.2fm]\n\n', r_min, r_max);

%% 3. 运行 GP 联合优化 (继承一维框架)
% ---------------------------------------------------------
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

% (1) 运行 NIGP (处理投影后的一维高斯)
theta0_nigp = [4.0, 2.0, 0.5, 5.0]; 
lb_nigp = [3, 0, 0.01, 0.1];
ub_nigp = [8, 3, 2.0,  20.0];
obj_nigp = @(theta) gp_nmll_nigp(theta, x_old, y_old_noisy, x1_exact, y1_exact, mu_r, yfar_exact, sigma_old, sigma_new, sigma_r);
opt_nigp = fmincon(obj_nigp, theta0_nigp, [], [], [], [], lb_nigp, ub_nigp, [], options);
A_nigp = opt_nigp(1); B_nigp = opt_nigp(2);

% (2) 运行 GP-LVM (处理投影后的一维区间)
theta0_lvm = [4.0, 2.0, (r_min+r_max)/2, 0.5, 5.0]; 
lb_lvm = [3, 0, r_min, 0.01, 0.1];
ub_lvm = [8, 3, r_max, 2.0,  20.0];
obj_lvm = @(theta) gp_nmll_lvm(theta, x_old, y_old_noisy, x1_exact, y1_exact, yfar_exact, sigma_old, sigma_new);
opt_lvm = fmincon(obj_lvm, theta0_lvm, [], [], [], [], lb_lvm, ub_lvm, [], options);
A_lvm = opt_lvm(1); B_lvm = opt_lvm(2); r_lvm = opt_lvm(3);

fprintf('=== 最终拟合结果 ===\n');
fprintf('真实参数: A = %.4f, B = %.4f\n', A_true, B_true);
fprintf('NIGP (高斯): A = %.4f, B = %.4f\n', A_nigp, B_nigp);
fprintf('GP-LVM(均匀): A = %.4f, B = %.4f, 估计真实距离 = %.4fm\n', A_lvm, B_lvm, r_lvm);

%% 4. 可视化
figure('Position', [100, 100, 1000, 450]);

% 子图 1: 2D 几何空间分布图
subplot(1, 2, 1); hold on; grid on;
% 画出二维矩形边界
patch([X_min X_max X_max X_min], [Y_min Y_min Y_max Y_max], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'm', 'DisplayName', '2D 均匀矩形');
% 画出二维高斯中心
plot(mu_2d(1), mu_2d(2), 'b*', 'MarkerSize', 8, 'DisplayName', '2D 高斯中心');
% 画射线
plot([0, 40*cos(theta_ray)], [0, 40*sin(theta_ray)], 'k--', 'LineWidth', 1.5, 'DisplayName', '1D 测量射线');
% 画真实位置
plot(pos_true_2d(1), pos_true_2d(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', '真实位置');
% 画均匀截线
plot([r_min*cos(theta_ray), r_max*cos(theta_ray)], [r_min*sin(theta_ray), r_max*sin(theta_ray)], 'm-', 'LineWidth', 4, 'DisplayName', '截获 1D 区间');
% 画高斯截点
plot(mu_r*cos(theta_ray), mu_r*sin(theta_ray), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '1D 高斯等效均值');

axis equal; xlim([0, 40]); ylim([0, 25]);
xlabel('X 坐标 (m)'); ylabel('Y 坐标 (m)'); title('二维不确定性几何切片');
legend('Location', 'best');

% 子图 2: 1D 信号强度拟合图
subplot(1, 2, 2); hold on; grid on;
x_plot = linspace(0.5, 45, 200)';
plot(x_plot, A_true + B_true*log(x_plot), 'k-', 'LineWidth', 2, 'DisplayName', '真实趋势');
scatter(x_old, y_old_noisy, 10, [0.7 0.7 0.7], 'filled', 'DisplayName', '带空间相关性历史数据');
plot(x_plot, A_nigp + B_nigp*log(x_plot), 'b--', 'LineWidth', 2, 'DisplayName', 'NIGP (2D->1D 高斯)');
plot(x_plot, A_lvm + B_lvm*log(x_plot), 'm-.', 'LineWidth', 2, 'DisplayName', 'GP-LVM (2D->1D 均匀)');
plot(r_lvm, yfar_exact, 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'GP-LVM 最终定位');

xlabel('射线一维距离 r (m)'); ylabel('信号强度 y'); title('GP 回归模型拟合效果');
legend('Location', 'best');

%% 辅助函数 (同之前一致)
function nmll = gp_nmll_nigp(theta, X_old, Y_old, x1, y1, x_far, y_far, sig_old, sig_new, sig_x)
    A = theta(1); B = theta(2); sig_f = theta(3); l = theta(4);
    X_all = [X_old; x1; x_far]; Y_all = [Y_old; y1; y_far]; N = length(Y_all);
    mu = A + B * log(X_all); res = Y_all - mu;
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    var_far_inflated = sig_new^2 + (B / x_far)^2 * sig_x^2;
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; var_far_inflated]);
    C = K + Noise;
    try
        L_chol = chol(C, 'lower'); alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch, nmll = Inf; end
end

function nmll = gp_nmll_lvm(theta, X_old, Y_old, x1, y1, y_far, sig_old, sig_new)
    A = theta(1); B = theta(2); x_far = theta(3); sig_f = theta(4); l = theta(5);
    X_all = [X_old; x1; x_far]; Y_all = [Y_old; y1; y_far]; N = length(Y_all);
    mu = A + B * log(X_all); res = Y_all - mu;
    K = sig_f^2 * exp(-0.5 * (pdist2(X_all, X_all) / l).^2);
    Noise = diag([sig_old^2 * ones(length(X_old), 1); sig_new^2; sig_new^2]);
    C = K + Noise;
    try
        L_chol = chol(C, 'lower'); alpha = L_chol' \ (L_chol \ res);
        nmll = 0.5 * res' * alpha + sum(log(diag(L_chol))) + 0.5 * N * log(2*pi);
    catch, nmll = Inf; end
end