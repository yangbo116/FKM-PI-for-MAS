% dimension 
% n = 8;
n = 12;
% n = 16;
% n = 20;
% n = 24;

% data
beta_1 = 1/(2*n + 2);
t1 = 6;
t2 = -1 - beta_1;
t3 = -1 + beta_1;

M2 = diag(ones(1, n-1), -1);
M1 = diag(ones(1, n),    0);
M3 = diag(ones(1, n-1),  1);

T1 = t2*M2 + t1*M1 + t3*M3;
T2 = T1 - t1*M1;
T3 = T2;

% construct linear system A
A1 = kron(T1, eye(n));
A  = kron(A1, eye(n));
A2 = kron(eye(n), T2);
B  = kron(A2, eye(n));
A3 = kron(eye(n), eye(n));
C  = kron(A3, T3);
A  = A + B + C;

n  = power(n,3);
b  = A*ones(n, 1);

% Parameter setting
tol    = 1e-6;
maxit  = 1000;
x_init = zeros(n, 1);

% HSS 
% n=8,  
% alphaADI=2.0521; omegaADI=0.0;
% n=12, 
alphaADI=1.4359;omegaADI=0.0;
% n=16,
% alphaADI=1.1025;omegaADI=0.0;
% n=20, 
% alphaADI=0.8943;omegaADI=0.0;
% n=24, 
% alphaADI=0.7520;omegaADI=0.0;

ADI_tim = tic;
[k_ADI, ~, Res_HSS] = GADI_HS(A, b, alphaADI,omegaADI, tol, maxit, x_init);
ADI_end = toc(ADI_tim);

% GADI_HS 
% n=8,  
% alphaGADI=0.6208;omegaGADI=1.0;
% n=12, 
alphaGADI=0.4468;omegaGADI=1.0;
% n=16, 
% alphaGADI=0.3465;omegaGADI=1.0;
% n=20, 
% alphaGADI=0.2823;omegaGADI=1.0;
% n=24,
% alphaGADI=0.2380;omegaGADI=1.0;

GADI_tim = tic;
[k_GADI, ~, Res_GAD_HS] = GADI_HS(A, b, alphaGADI, omegaGADI, tol, maxit, x_init);
GADI_end = toc(GADI_tim);

% AGADI_HS = GAD_HS + halpern, omga\in[0,2) convergence
% omegaAGADI_HS = 0: AGADI_HS--->PRS + halpern
% n=8,  
% alphaAGADI=0.50; omegaAGADI=0.62;
% n=12, 
alphaAGADI=0.29; omegaAGADI=0.64;
% n=16, 
% alphaAGADI=0.20; omegaAGADI=0.64;
% n=20, 
% alphaAGADI=0.17; omegaAGADI=0.54;
% n=24, 
% alphaAGADI=0.13; omegaAGADI=0.50;

AGADI_tim = tic;
[k_AGADI, ~, Res_AGADI_HS] = AGADI_HS(A, b, alphaAGADI, omegaAGADI, tol, maxit, x_init);
AGADI_end = toc(AGADI_tim);

% plot HSS, GADI_HS, AGADI_HS
% figure;
% HSS
semilogy(Res_HSS, '-+','LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', ...
    'ADI(\alpha^{\ast}='+ string(alphaADI)+')');
hold on;

% GADI_HS
semilogy(Res_GAD_HS, '-+','LineWidth',1.5, 'MarkerSize', 5, 'DisplayName', ...
     'GADI(\alpha^{\ast}='+ string(alphaGADI)+', \omega^{\ast}='+ string(omegaGADI)+')');

% AGADI_HS
semilogy(Res_AGADI_HS, '-+','LineWidth',1.5, 'MarkerSize', 5, 'DisplayName', ...
     'FKM\_ADI(\alpha^{\ast}='+ string(alphaAGADI)+', \omega^{\ast}='+ string(omegaAGADI)+')');

% OC_AGADI
% semilogy(OC_R_res, '-+','LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', ...
    % 'OC\_AGADI\_HS(a='+ string(OCalpAGADI)+', w='+ string(OComegaAGADI)+')');

xlabel('Iteration number');
ylabel('Residual=$\frac{\|r^{k}\|_{2}}{\|r^{0}\|_{2}}$','Interpreter', 'latex','FontSize', 16);
title('n = '+ string(n),'Interpreter', 'latex', 'FontSize', 13);
grid on;
legend('Location','best');
axis tight;
hold off;

fprintf('ADI:   CPU time: %.4f,iteration number: %d \n',  ADI_end, k_ADI);
fprintf('GADI:  CPU time: %.4f,iteration number: %d \n', GADI_end, k_GADI);
fprintf('AGDAI: CPU time: %.4f,iteration number: %d \n', AGADI_end, k_AGADI);