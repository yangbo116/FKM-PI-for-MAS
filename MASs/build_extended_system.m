function [A_bar, B_bar, C_bar, Q_bar] = build_extended_system(A, B, C, F, T, Q)
    ni = size(A, 1);
    p = size(F, 1);
    A_bar = blkdiag(A, F);
    B_bar = [B; zeros(p, size(B, 2))];
    C_bar = blkdiag(C, T);
    Q_bar = [C' * Q * C, -C' * Q * T; -T' * Q * C, T' * Q * T];
end