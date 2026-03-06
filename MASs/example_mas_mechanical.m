clear; clc;
N = 6;

km=1; ju=0.2; bm=0.2; lu=0.5;

A{1} = [0, 1; -km/ju, -bm/ju];  B{1} = [0; 1/lu];  C{1} = [1, 0];
A{2} = [0, 1; -km/ju, -bm/ju];  B{2} = [0; 1/lu];  C{2} = [1, 0];
A{3} = [0, 1; -km/ju, -bm/ju];  B{3} = [0; 1/lu];  C{3} = [1, 0];
A{4} = [0, 1; -km/ju, -bm/ju];  B{4} = [0; 1/lu];  C{4} = [1, 0];
A{5} = [0, 1, 0; 0, 0, 1; 0, 0, -0.7];  B{5} = [6; 0; 1];  C{5} = [1, -1, 1];
A{6} = [0, 1, 0; 0, 0, 1; 0, 0, -0.7];  B{6} = [6; 0; 1];  C{6} = [1, -1, 1];

F = [0, 1; -1, -0.5]; 
T = [1, 0];
Q = 2*eye(1);
R = eye(1);
epsilon = 1e-8;
p = size(F, 1);



X_0 = [];
for i = 1:N
    ni = size(A{i},1);                
    agent_initial_state = 0.5*i * 1e2 * ones(ni, 1); 
    X_0 = [X_0; agent_initial_state];  
end



X_0 = [X_0; 1e2*ones(p*N, 1); 1e2*ones(p, 1)];


A_adj = [0 0 0 0 0 0;
         0 0 0 0 0 1;
         1 0 0 1 0 1;
         0 0 1 0 0 1;
         1 0 0 0 0 1;
         0 0 1 0 1 0];
G = diag([1, 0, 0, 0, 1, 1]);
c = 0.5;
dt = 1;

K_opt  = cell(1, N);
K_iter = cell(1, N);
P_opt  = cell(1, N);
P_iter = cell(1, N);
P_norms = cell(1, N);
Y_hist  = cell(1, N);
E_hist  = cell(1, N);

K_MAX  = 20;
K_time = 50;
k_values = 0:K_time;


for i = 1:N
    
    [A_bar, B_bar, C_bar, Q_bar] = build_extended_system(A{i}, B{i}, C{i}, F, T, Q);

    [P_opt{i}, K_opt{i}] = solve_riccati_optimal(A_bar, B_bar, Q_bar, R);

    [K_iter{i}, P_iter{i}, P_norms{i}, K_norms{i}] = solve_iterative_gain(A_bar, B_bar, Q_bar, R, epsilon, P_opt{i}, K_opt{i},K_MAX);

    Y_hist{i} = [];
    E_hist{i} = [];
    U_hist{i} = [];

    for k = k_values
        [X_k, Y_k, u_k, e_k, R0_k] = compute_states_error( ...
            A_bar, B_bar, C_bar, K_iter{i}, X_0, k, A_adj, G, c, F, T, dt, i, A);

        R0_k_agents = reshape(R0_k(1:N*p), p, N);
        y0_i_k = T * R0_k_agents(:, i);

        e_k = norm(Y_k(1) - y0_i_k, 2);

        Y_hist{i} = [Y_hist{i}, Y_k];
        E_hist{i} = [E_hist{i}, e_k];
        U_hist{i} = [U_hist{i}, u_k];
    end
end


agent_leg = arrayfun(@(i) sprintf('Agent%d', i), 1:N, 'UniformOutput', false);
y_leg     = arrayfun(@(i) sprintf('y_%d', i),     1:N, 'UniformOutput', false);
e_leg     = arrayfun(@(i) sprintf('e_%d', i),     1:N, 'UniformOutput', false);
u_leg     = arrayfun(@(i) sprintf('u_%d', i),     1:N, 'UniformOutput', false);
clr = lines(N);


figure; hold on;
for i = 1:N
    plot(0:(length(K_norms{i})-1), K_norms{i}, 'LineWidth', 2, 'Color', clr(i,:));
end
xlabel('Iteration number','FontSize',16);
ylabel('$$\|K_i(k)-K_{opt}\|_2$$','Interpreter','latex','FontSize',16);
legend(agent_leg{:},'FontSize',16,'Location','best');
set(gca,'LineWidth',2,'FontSize',16,'Box','on'); grid on; hold off;


figure; hold on;
for i = 1:N
    plot(0:(length(P_norms{i})-1), P_norms{i}, 'LineWidth', 2, 'Color', clr(i,:));
end
xlabel('Iteration number','FontSize',16);
ylabel('$$\|P_i(k)-P_{opt}\|_2$$','Interpreter','latex','FontSize',16);
legend(agent_leg{:},'FontSize',16,'Location','best');
set(gca,'LineWidth',2,'FontSize',16,'Box','on'); grid on; hold off;


figure; hold on;
for i = 1:N
    
    plot(k_values, Y_hist{i}(1,:), 'LineWidth', 2, 'Color', clr(i,:));
end

plot(k_values, Y_hist{1}(end,:), 'k--', 'LineWidth', 3);

xlabel('Time(s)','FontSize',16);
ylabel('Output trajectories','FontSize',16);
legend([y_leg, {'y_0'}], 'FontSize',16,'Location','best');
set(gca,'LineWidth',2,'FontSize',16,'Box','on'); grid on; hold off;


figure; hold on;
for i = 1:N
    plot(k_values, E_hist{i}, 'LineWidth', 2, 'Color', clr(i,:));
end
xlabel('Time(s)','FontSize',16);
ylabel('$$\|e_i(k)\|_2$$','Interpreter','latex','FontSize',16);
legend(e_leg{:},'FontSize',16,'Location','best');
set(gca,'LineWidth',2,'FontSize',16,'Box','on'); grid on; hold off;


figure; hold on;
for i = 1:N
    plot(k_values, U_hist{i}, 'LineWidth', 2, 'Color', clr(i,:));
end
xlabel('Time(s)','FontSize',16);
ylabel('$$u_i(k)$$','Interpreter','latex','FontSize',16);
legend(u_leg{:},'FontSize',16,'Location','best');
set(gca,'LineWidth',2,'FontSize',16,'Box','on'); grid on; hold off;