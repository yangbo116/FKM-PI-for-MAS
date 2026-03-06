% Data sources: On inexact alternating direction implicit iteration for
% continuous Sylvester equations, example 1

n = 256;
r = 1;

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
tol   = 1e-6;
maxit = 1000;
X_init= zeros(n,n);

GADIomge = 0.1;
alp1 = 0.54; beta1=0.54;
GADI_tim = tic;
[k_GADI, ~, R_GADI] = GADI_AB(A, B, C, alp1, beta1, GADIomge, tol, maxit, X_init);
GADI_end = toc(GADI_tim);

AGADIomge = 0.01;
alp = 0.43; beta=0.43;
AGADI_tim = tic;
[k_AGADI, ~, R_AGADI] = AGADI_AB(A, B, C, alp, beta, AGADIomge, tol, maxit, X_init);
AGADI_end = toc(AGADI_tim);

% ADI
figure;
semilogy(R_GADI, '-+','LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', ...
    'GADI(a='+ string(alp1)+', w='+ string(GADIomge)+')');
hold on;

% FKMADI 
semilogy(R_AGADI, '-+', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', ...
    'AGADI(a='+ string(alp)+', w='+ string(AGADIomge)+')');

% Plot 
xlabel('Iteration number');
ylabel('Residual=$\frac{\|R^{k}\|_{F}}{\|C\|_{F}}$','Interpreter', 'latex','FontSize', 16);
title('Performance comparison for n='+string(n) +', $r=$'+string(r) ...
    , 'Interpreter', 'latex', 'FontSize', 13);
grid on;
legend('Location','best');
axis tight;
hold off;

fprintf('GDAI:  CPU time: %.4f,iteration number: %d \n',  GADI_end, k_GADI);
fprintf('AGDAI: CPU time: %.4f,iteration number: %d \n', AGADI_end, k_AGADI);