% =========================================================
% Constrained WLS Fitting for A + B*log(x)
% 带约束的加权最小二乘法曲线拟合
% =========================================================
clear; clc; close all;

%% 1. 参数设置 (手动设置区域)
% ---------------------------------------------------------
% 真实值 (用于生成新版本的“准”观测点，也是我们拟合要逼近的终极目标)
A_true = 5.5; 
B_true = 1.2;

% 老版本曲线参数 (用于生成历史的、带噪声的数据集基准)
A_old = 4.0;
B_old = 2.0;
noise_std = 0.5; % 历史数据的噪声标准差

% 约束范围设定 (以老版本为基准的上下波动幅度)
delta_A = 3.0; % A的浮动范围
delta_B = 1.5; % B的浮动范围

% 准的训练点 x 位置
x0 = 5; 
x_new = [1; x0]; % x=1 和 x=x0

% 权重分配
W_new = 100; % 新观测点极其准确，赋予极高权重
W_old = 1;   % 老数据集主要起基准平滑作用，权重较低

%% 2. 生成数据集
% ---------------------------------------------------------
% 2.1 生成老版本带噪声的曲线库 (x 范围从 0.5 到 10)
N_old = 150;
x_old = linspace(0.5, 10, N_old)';
y_old_clean = A_old + B_old * log(x_old);
y_old_noisy = y_old_clean + noise_std * randn(N_old, 1);

% 2.2 生成最新的“准”观测点 (基于真实值，无噪声)
y_new = A_true + B_true * log(x_new);

%% 3. 构建并求解带约束的加权最小二乘 (WLS)
% ---------------------------------------------------------
% 目标函数形式: || C * theta - d ||^2，其中 theta = [A; B]
% 注意：为了匹配 lsqlin 的目标函数形式，矩阵需乘以 sqrt(W)

% 构建老数据的设计矩阵和目标向量
C_old = sqrt(W_old) * [ones(N_old, 1), log(x_old)];
d_old = sqrt(W_old) * y_old_noisy;

% 构建新数据的设计矩阵和目标向量
C_new = sqrt(W_new) * [ones(length(x_new), 1), log(x_new)];
d_new = sqrt(W_new) * y_new;

% 合并矩阵
C = [C_new; C_old];
d = [d_new; d_old];

% 设定上下界 (Hard Constraints: [lb, ub])
lb = [A_old - delta_A; B_old - delta_B];
ub = [A_old + delta_A; B_old + delta_B];

% 调用 lsqlin 求解优化问题
options = optimoptions('lsqlin', 'Display', 'off');
theta_opt = lsqlin(C, d, [], [], [], [], lb, ub, [], options);

A_opt = theta_opt(1);
B_opt = theta_opt(2);

% 打印结果
fprintf('=== 拟合结果 ===\n');
fprintf('真实参数: A = %.4f, B = %.4f\n', A_true, B_true);
fprintf('拟合参数: A = %.4f, B = %.4f\n', A_opt, B_opt);
fprintf('优化误差: A误差 = %.4f, B误差 = %.4f\n', abs(A_true-A_opt), abs(B_true-B_opt));

%% 4. 可视化对比
% ---------------------------------------------------------
x_plot = linspace(0.5, 10, 500)';
y_true_plot = A_true + B_true * log(x_plot);
y_old_plot  = A_old + B_old * log(x_plot);
y_opt_plot  = A_opt + B_opt * log(x_plot);

figure('Position', [100, 100, 800, 600]);
hold on; grid on;

% 绘制老版本的噪声数据集
scatter(x_old, y_old_noisy, 20, [0.7 0.7 0.7], 'filled', 'DisplayName', '旧版本噪声数据');

% 绘制旧版本基准曲线
plot(x_plot, y_old_plot, 'g--', 'LineWidth', 1.5, 'DisplayName', '旧版本基准曲线');

% 绘制真实的新版本曲线
plot(x_plot, y_true_plot, 'k-', 'LineWidth', 2, 'DisplayName', '新版本真实曲线 (Ground Truth)');

% 绘制WLS拟合出的曲线
plot(x_plot, y_opt_plot, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Constrained WLS 拟合曲线');

% 绘制两个“准”的训练点
plot(x_new, y_new, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', '高权重准点');

% 图形装饰
xlabel('x');
ylabel('y');
title('带约束的加权最小二乘法 (Constrained WLS) 曲线拟合');
legend('Location', 'best');
set(gca, 'FontSize', 12);