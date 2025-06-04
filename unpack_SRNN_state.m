function [a, b, u_d] = unpack_SRNN_state(X, n, n_a, n_b)
    % Unpack SRNN state vector/matrix X into named variables.
    % If X is N_sys_eqs x 1 (single time point, e.g., from within SRNN.m):
    %   a will be n x n_a
    %   b will be n x n_b
    %   u_d will be n x 1
    % If X is nt x N_sys_eqs (multiple time points, e.g., output from ode45):
    %   a will be n x n_a x nt
    %   b will be n x n_b x nt
    %   u_d will be n x nt

    is_single_time_point = iscolumn(X);

    if is_single_time_point
        nt = 1; % Conceptually, for internal calculations
    else
        nt = size(X, 1); % Number of time points from the matrix X
    end

    idx_a_end = n * n_a;

    % ---------- SFA (a–states) ----------
    if n_a > 0
        if is_single_time_point
            % X is (N_sys_eqs) x 1. X(1:idx_a_end) is (n*n_a) x 1.
            a = reshape(X(1:idx_a_end), n, n_a);      % n x n_a
        else
            % X is nt x N_sys_eqs. X_a_block is nt x (n*n_a).
            X_a_block = X(:, 1:idx_a_end);
            X_a_block_T = X_a_block'; % (n*n_a) x nt
            a = reshape(X_a_block_T, n, n_a, nt); % n x n_a x nt
        end
    else
        a = []; % no SFA states
    end

    % ---------- STD (b–states) ----------
    idx_b_start = idx_a_end + 1;
    idx_b_end   = idx_a_end + n * n_b;

    if n_b > 0
        if is_single_time_point
            % X_b_component is (n*n_b) x 1
            b = reshape(X(idx_b_start:idx_b_end), n, n_b); % n x n_b
        else
            % X_b_block is nt x (n*n_b)
            X_b_block = X(:, idx_b_start:idx_b_end);
            X_b_block_T = X_b_block'; % (n*n_b) x nt
            b = reshape(X_b_block_T, n, n_b, nt); % n x n_b x nt
        end
    else
        b = []; % no STD states
    end

    % ---------- Dendrite (u_d) ----------
    idx_ud_start = idx_b_end + 1;
    % If single time point, idx_ud_end is total length of X (N_sys_eqs).
    % If multiple time points, idx_ud_end is number of columns in X (N_sys_eqs).
    idx_ud_end = idx_b_end + n;

    if is_single_time_point
        % X_ud_component is n x 1
        u_d = X(idx_ud_start:idx_ud_end); % n x 1
    else
        % X_ud_block is nt x n
        X_ud_block = X(:, idx_ud_start:idx_ud_end);
        u_d = X_ud_block'; % n x nt
    end
end