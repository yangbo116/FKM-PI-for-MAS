function [k, x, R_res] = AGADI_HS(A, b, alpha, omega, tol, maxit, x_init)

% x = dA\b solving Ax = b
% x = b/dA solving xA = b

[m, n] = size(A);

% contruct matrix H and S
H = (A + A')/2;
S = (A - A')/2;

% dH cach the inverse of alpha*eye(m) + H
% dS cach the inverse of alpha*eye(n) + H

dH = decomposition(alpha*eye(m)+ H);
dS = decomposition(alpha*eye(n)+ S);

% 1、(alpha*I+H)x = rhs1，
% 2、(alpha*I+S)x = rhs2

rhs_1 = alpha*eye(m) - S;
rhs_2 = S - (1 - omega)*alpha*eye(n);

R_res = zeros(maxit,1);
R_res_init = norm(A*x_init - b);
x_old = x_init;
l = 1;
for k = 1: maxit
    
    x = dH\(rhs_1*x_old + b);
    x_hat = dS\(rhs_2*x_old + (2-omega)*alpha*x);
    
    % Halpern iteration update with fixed step restart
    x = (1/(l + 30))*x_init + ((l + 29)/(l + 30))*x_hat;

    if mod(k, 2)==0
        l = 1;
        x_init = x_hat;
        x_old  = x_hat;
    else
        l = l + 1;
        x_old = x;
    end

    res_norm = norm(A*x-b)/R_res_init;
    R_res(k) = res_norm;
    if res_norm < tol
        R_res = R_res(1:k);
        break
    end
end

end

