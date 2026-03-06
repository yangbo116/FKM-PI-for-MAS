function [K, P, P_norms, K_norms] = solve_iterative_gain(A_bar, B_bar, Q_bar, R, epsilon, P_opt, K_opt, max_iter)
    mi = size(B_bar, 2);
    n = size(A_bar, 1);
    % K = 1e2 * ones(mi, n);
    K = 1e2 * rand(mi, n);
    P = eye(n);
    P_norms = [];
    K_norms = [];

    for iter = 1:max_iter
        P_next = sylvester((A_bar - B_bar * K)',(A_bar - B_bar * K),-(Q_bar + K' * R * K));
        K_next = inv(R) * B_bar' * P_next;
        P_norms = [P_norms, norm(P_next - P_opt, 2)];
        K_norms = [K_norms, norm(K_next - K_opt, 2)];

        if norm(K_next - K, 'fro') <= epsilon
            break;
        end
        
        P = P_next;
        K = K_next;
    end
end
