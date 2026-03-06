function [X_k, Y_k, u_k, e_k, R0_k] = compute_states_error( ...
    A_bar, B_bar, C_bar, K, X_0, k, A_adj, G, c, F, T, dt, idx_agent, A_cell)

    N = size(A_adj, 1);
    p = size(F, 1);

    n_list = cellfun(@(Ai) size(Ai,1), A_cell);     
    idx0   = cumsum([1, n_list(1:end-1)]);           
    ni     = n_list(idx_agent);
    x0_i = X_0(idx0(idx_agent) : idx0(idx_agent)+ni-1);
    R0_initial = X_0(end - (N+1)*p + 1 : end);
    R0_agents0 = reshape(R0_initial(1:N*p), p, N);
    eta0_i = R0_agents0(:, idx_agent);
    function dRdt = odefun(~, R)
        R0_agents = reshape(R(1:N*p), p, N);
        R0_global = R(N*p + 1 : end);

        coupling = zeros(p, N);
        for ii = 1:N
            sum_neighbors = zeros(p, 1);
            for jj = 1:N
                sum_neighbors = sum_neighbors + A_adj(ii,jj) * (R0_agents(:,jj) - R0_agents(:,ii));
            end
            coupling(:,ii) = c * (sum_neighbors - G(ii,ii) * (R0_global - R0_agents(:,ii)));
        end

        dR0_agents = F * R0_agents + coupling;
        dR0_global = F * R0_global;

        dRdt = [dR0_agents(:); dR0_global(:)];
    end
    if k == 0
        R0_k = R0_initial;

        X_k = [x0_i; eta0_i];  
        assert(size(X_k,1) == size(C_bar,2), "Dimension mismatch: X_k vs C_bar.");
        u_k = -K * X_k;

        Y_k = C_bar * X_k;    
        e_k = Y_k(1) - (T * eta0_i); 
        return;
    end

    tspan = [0, k*dt];
    [~, R] = ode45(@odefun, tspan, R0_initial);
    R0_k = R(end, :)';

    
    X_0_original = [x0_i; eta0_i];    
    assert(size(X_0_original,1) == size(A_bar,1), "Dimension mismatch: X_0_original vs A_bar.");

    Ak_bar = A_bar - B_bar * K;
    X_k = expm(Ak_bar * (k*dt)) * X_0_original;
    u_k = -K * X_k;
    Y_k = C_bar * X_k;

    R0_k_agents = reshape(R0_k(1:N*p), p, N);
    y0_i_k = T * R0_k_agents(:, idx_agent);

    y_i_k = Y_k(1);        
    e_k = y_i_k - y0_i_k;    
end

