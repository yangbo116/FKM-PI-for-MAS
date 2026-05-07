function [k, X, R_res] = GADI_AB(A, B, C, hatalp, beta, hatomg, tol, maxit, X_init)

% x = dA\b solving Ax = b
% x = b/dA solving xA = b
% Solving the Sylvester equation AX+XB=C by ADI

[m, n] = size(C);
dA = decomposition(hatalp*eye(m) + A);
dB = decomposition(beta*eye(n) + B);

B_alpha = hatalp*eye(n) - B;
A_beta  = B - (1-hatomg)*hatalp*eye(n);

X_old = X_init;

norm_C= norm(C,'fro');
R_res = zeros(maxit,1);

for k = 1:maxit

    X_half= dA\(X_old*B_alpha + C);
    X = (X_old*A_beta + (2-hatomg)*hatalp*X_half)/dB;

    R = C - A*X - X*B;
    res_norm = norm(R,'fro')/norm_C;
    R_res(k) = res_norm;
    if res_norm < tol
        R_res = R_res(1:k);
        break
    else
        X_old = X;
    end
end

end
