% Data sources: On inexact alternating direction implicit iteration for
% continuous Sylvester equations, example 1

n = 256;
r = 0.01;

% Prepare A, B, C
rng(142);
M1= diag(2*ones(1, n),0);
M2= diag(ones(1, n-1),-1);
M3= diag(ones(1, n-1),1);
M = (-1)*M2 + M1 + (-1)*M3;
N = 0.5*M2 + zeros(n) + (-0.5)*M3;
A = M+2*r*N+(100/(n+1)^2)*eye(n);
B = A;
X = randn(n, n);
C = A*X+X*B;

% Parameter setting
tol = 1e-6;
maxit = 1000;
X_init= zeros(n,n);

% ==========================================================
% 第一阶段：粗网格并行搜索 (Coarse Grid Search)
% ==========================================================
fprintf('--- 开始第一阶段：粗网格快速搜索 ---\n');

% 粗网格步长 0.1
alpha_coarse = 0.01 : 0.1 : 3.0; 
omega_coarse = 0.00 : 0.1 : 1.99;

num_a_c = length(alpha_coarse);
num_o_c = length(omega_coarse);

Iter_Matrix_C = zeros(num_a_c, num_o_c);
FinalRes_Matrix_C = zeros(num_a_c, num_o_c);

% 开启并行池 (如果未开启的话)
if isempty(gcp('nocreate')), parpool; end

% 使用 parfor 进行多核并行计算
parfor i = 1:num_a_c
    alpha = alpha_coarse(i);
    beta = alpha; 
    
    % parfor 内部的临时行变量
    temp_iter = zeros(1, num_o_c);
    temp_res = zeros(1, num_o_c);
    
    for j = 1:num_o_c
        omega = omega_coarse(j);
        [k, ~, res] = AGADI_AB(A, B, C, alpha, beta, omega, tol, maxit, X_init);
        temp_iter(j) = k;
        temp_res(j) = res(end);
    end
    Iter_Matrix_C(i, :) = temp_iter;
    FinalRes_Matrix_C(i, :) = temp_res;
end

% 提取粗搜索最优解（包含处理平局的逻辑）
min_k_C = min(Iter_Matrix_C(:));
idx_C = find(Iter_Matrix_C == min_k_C);
[~, best_tie_C] = min(FinalRes_Matrix_C(idx_C));
[best_row_C, best_col_C] = ind2sub([num_a_c, num_o_c], idx_C(best_tie_C));

opt_alpha_coarse = alpha_coarse(best_row_C);
opt_omega_coarse = omega_coarse(best_col_C);

fprintf('粗搜完成！大致最优区域: alpha = %.2f, omega = %.2f (步数: %d)\n\n', ...
    opt_alpha_coarse, opt_omega_coarse, min_k_C);

% ==========================================================
% 第二阶段：局部细网格并行搜索 (Fine Grid Search)
% ==========================================================
fprintf('--- 开始第二阶段：局部细网格精确搜索 ---\n');

% 在粗搜最优点附近 +/- 0.15 的范围内，使用 0.01 步长细搜
% 注意：边界不能超出原始定义的 [0.01, 3.0] 和 [0.00, 1.99]
alpha_fine = max(0.01, opt_alpha_coarse - 0.15) : 0.01 : min(3.0, opt_alpha_coarse + 0.15);
omega_fine = max(0.00, opt_omega_coarse - 0.15) : 0.01 : min(1.99, opt_omega_coarse + 0.15);

num_a_f = length(alpha_fine);
num_o_f = length(omega_fine);

Iter_Matrix_F = zeros(num_a_f, num_o_f);
FinalRes_Matrix_F = zeros(num_a_f, num_o_f);

parfor i = 1:num_a_f
    alpha = alpha_fine(i);
    beta = alpha; 
    
    temp_iter = zeros(1, num_o_f);
    temp_res = zeros(1, num_o_f);
    
    for j = 1:num_o_f
        omega = omega_fine(j);
        [k, ~, res] = AGADI_AB(A, B, C, alpha, beta, omega, tol, maxit, X_init);
        temp_iter(j) = k;
        temp_res(j) = res(end);
    end
    Iter_Matrix_F(i, :) = temp_iter;
    FinalRes_Matrix_F(i, :) = temp_res;
end

% 提取细搜索绝对最优解
min_k_F = min(Iter_Matrix_F(:));
idx_F = find(Iter_Matrix_F == min_k_F);
[min_final_res, best_tie_F] = min(FinalRes_Matrix_F(idx_F));
[best_row_F, best_col_F] = ind2sub([num_a_f, num_o_f], idx_F(best_tie_F));

opt_alpha = alpha_fine(best_row_F);
opt_omega = omega_fine(best_col_F);

fprintf('============== 最终寻优结果 ==============\n');
fprintf('最优 alpha: %.4f\n', opt_alpha);
fprintf('最优 omega: %.4f\n', opt_omega);
fprintf('最少迭代步数: %d\n', min_k_F);
fprintf('最优参数下的最终残差: %e\n', min_final_res);
fprintf('==========================================\n');

% ==========================================================
% 第三阶段：可视化 (数据画图)
% ==========================================================
fprintf('--- 正在生成参数空间寻优可视化图表 ---\n');

% 创建一个宽幅图窗
figure('Name', '粗细网格两阶段参数寻优可视化', 'Position', [100, 100, 1100, 450]);
colormap('parula'); % 设置统一的色彩映射

% ----------------------------------------------------------
% 子图 1：粗网格全局搜索等高线图
% ----------------------------------------------------------
subplot(1, 2, 1);
[Omega_Grid_C, Alpha_Grid_C] = meshgrid(omega_coarse, alpha_coarse);

% 绘制粗网格等高线（填充颜色）
contourf(Omega_Grid_C, Alpha_Grid_C, Iter_Matrix_C, 25, 'LineColor', 'none');
colorbar; hold on;

% 标出粗搜找到的最优点 (蓝点)
plot(opt_omega_coarse, opt_alpha_coarse, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% 绘制细搜区域的边界框 (红色虚线框)
w_min = omega_fine(1); w_max = omega_fine(end);
a_min = alpha_fine(1); a_max = alpha_fine(end);
rect_pos = [w_min, a_min, w_max - w_min, a_max - a_min];
rectangle('Position', rect_pos, 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 2);

% 设置坐标轴和标签
xlabel('\omega');
ylabel('\alpha');
title('Stage one：global coarse search');
legend('Iteration contours', 'Coarse Optimum', 'Fine Zoom Region', 'location', 'northeast');
axis tight; 

% ----------------------------------------------------------
% 子图 2：细网格局部精确搜索等高线图
% ----------------------------------------------------------
subplot(1, 2, 2);
[Omega_Grid_F, Alpha_Grid_F] = meshgrid(omega_fine, alpha_fine);

% 绘制细网格等高线
contourf(Omega_Grid_F, Alpha_Grid_F, Iter_Matrix_F, 15, 'LineColor', 'none');
colorbar; hold on;

% 找出所有达到最小步数 k_F 的点 (平局点)
[tie_rows, tie_cols] = find(Iter_Matrix_F == min_k_F);
plot(omega_fine(tie_cols), alpha_fine(tie_rows), 'g.', 'MarkerSize', 10);

% 标出经过残差决胜后的最终绝对最优点 (红星)
plot(opt_omega, opt_alpha, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% 设置坐标轴和标签
xlabel('\omega');
ylabel('\alpha');
title('Stage two：local fine search');
legend('Iteration', 'Tied min iterations', 'Optimal residual', 'Location', 'northeast');
axis tight;

% 美化整体排版
sgtitle('FKM\_ADI: Two-Stage grid search optimization results analysis for n = '+string(n) +', r = '+string(r)',  'FontSize', 14, 'FontWeight', 'bold');