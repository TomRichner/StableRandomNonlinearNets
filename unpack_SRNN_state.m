function [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X, params)
    % Unpack SRNN state vector/matrix X into named variables for E and I neurons.
    %
    % If X is N_sys_eqs x 1 (single time point, e.g., from within SRNN.m):
    %   a_E will be n_E x n_a_E (if n_a_E > 0, else [])
    %   a_I will be n_I x n_a_I (if n_a_I > 0, else [])
    %   b_E will be n_E x n_b_E (if n_b_E > 0, else [])
    %   b_I will be n_I x n_b_I (if n_b_I > 0, else [])
    %   u_d will be n x 1
    %
    % If X is nt x N_sys_eqs (multiple time points, e.g., output from ode solver):
    %   a_E will be n_E x n_a_E x nt (if n_a_E > 0, else [])
    %   a_I will be n_I x n_a_I x nt (if n_a_I > 0, else [])
    %   b_E will be n_E x n_b_E x nt (if n_b_E > 0, else [])
    %   b_I will be n_I x n_b_I x nt (if n_b_I > 0, else [])
    %   u_d will be n x nt

    n_E   = params.n_E;
    n_I   = params.n_I;
    n_a_E = params.n_a_E;
    n_a_I = params.n_a_I;
    n_b_E = params.n_b_E;
    n_b_I = params.n_b_I;
    n     = params.n;

    is_single_time_point = iscolumn(X);

    if is_single_time_point
        nt = 1; 
    else
        nt = size(X, 1); % Number of time points
    end

    current_idx = 0;

    % --- SFA states for E neurons (a_E) ---
    len_a_E = n_E * n_a_E;
    if n_a_E > 0 && n_E > 0
        if is_single_time_point
            % X_aE_vec is (n_E*n_a_E) x 1
            a_E_vec = X(current_idx + (1:len_a_E));
            a_E = reshape(a_E_vec, n_E, n_a_E); % n_E x n_a_E
        else
            % X_aE_block is nt x (n_E*n_a_E)
            X_aE_block = X(:, current_idx + (1:len_a_E));
            X_aE_block_T = X_aE_block'; % (n_E*n_a_E) x nt
            a_E = reshape(X_aE_block_T, n_E, n_a_E, nt); % n_E x n_a_E x nt
        end
    else
        a_E = []; % No SFA states for E neurons or no E neurons
    end
    current_idx = current_idx + len_a_E;

    % --- SFA states for I neurons (a_I) ---
    len_a_I = n_I * n_a_I;
    if n_a_I > 0 && n_I > 0
        if is_single_time_point
            a_I_vec = X(current_idx + (1:len_a_I));
            a_I = reshape(a_I_vec, n_I, n_a_I); % n_I x n_a_I
        else
            X_aI_block = X(:, current_idx + (1:len_a_I));
            X_aI_block_T = X_aI_block';
            a_I = reshape(X_aI_block_T, n_I, n_a_I, nt); % n_I x n_a_I x nt
        end
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;

    % --- STD states for E neurons (b_E) ---
    len_b_E = n_E * n_b_E;
    if n_b_E > 0 && n_E > 0
        if is_single_time_point
            b_E_vec = X(current_idx + (1:len_b_E));
            b_E = reshape(b_E_vec, n_E, n_b_E); % n_E x n_b_E
        else
            X_bE_block = X(:, current_idx + (1:len_b_E));
            X_bE_block_T = X_bE_block';
            b_E = reshape(X_bE_block_T, n_E, n_b_E, nt); % n_E x n_b_E x nt
        end
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;

    % --- STD states for I neurons (b_I) ---
    len_b_I = n_I * n_b_I;
    if n_b_I > 0 && n_I > 0
        if is_single_time_point
            b_I_vec = X(current_idx + (1:len_b_I));
            b_I = reshape(b_I_vec, n_I, n_b_I); % n_I x n_b_I
        else
            X_bI_block = X(:, current_idx + (1:len_b_I));
            X_bI_block_T = X_bI_block';
            b_I = reshape(X_bI_block_T, n_I, n_b_I, nt); % n_I x n_b_I x nt
        end
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;

    % --- Dendrite states (u_d) ---
    % u_d is always n x 1 (single time point) or n x nt (multiple time points)
    if is_single_time_point
        u_d = X(current_idx + (1:n)); % n x 1
    else
        X_ud_block = X(:, current_idx + (1:n)); % nt x n
        u_d = X_ud_block'; % n x nt
    end
end