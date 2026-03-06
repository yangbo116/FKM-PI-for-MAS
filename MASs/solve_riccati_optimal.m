function [P_opt, K_opt] = solve_riccati_optimal(A_bar, B_bar, Q_bar, R)
    try
        P_opt = care(A_bar, B_bar, Q_bar, R);
        K_opt = inv(R)* B_bar' *P_opt;
    catch
        P_opt = [];
        K_opt = [];
    end
end
