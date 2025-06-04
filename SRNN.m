function [dX_dt] = SRNN(~, X, u_ex, params)

    %% load parameters
    n_a = params.n_a; % number of SFA timescales per neuron
    n_b = params.n_b; % number of STD timescales per neuron

    tau_a = params.tau_a; % s, 1 x n_a, time constants of SFA
    tau_b = params.tau_b; % s, 1 x n_b, time constants of STD

    tau_d = params.tau_d; % s, scalar

    n = params.n_neurons;
    M = params.M;
    
    c_SFA = params.c_SFA; % n_neurons x 1, 0 for I neurons
    F_STD = params.F_STD; % n_neurons x 1, 0 for I neurons
    tau_STD = params.tau_STD; % scalar

    %% index state variables.
    idx_a_end = n * n_a;

    % ---------- SFA (a–states) ----------
    if n_a > 0
        a = reshape(X(1:idx_a_end), n, n_a);      % n_neurons × n_a
    else
        a = [];                                   % no SFA states
    end

    % ---------- STD (b–states) ----------
    idx_b_start = idx_a_end + 1;
    idx_b_end   = idx_a_end + n * n_b;

    if n_b > 0
        b = reshape(X(idx_b_start:idx_b_end), n, n_b); % n_neurons × n_b
    else
        b = [];                                   % no STD states
    end

    % ---------- Dendrite (u_d) ----------
    idx_ud_start = idx_b_end + 1;
    idx_ud_end   = idx_b_end + n;
    u_d = X(idx_ud_start:idx_ud_end);

    %% make dependent variables from state variables and parameters
    if n_a > 0
        r = relu(u_d - c_SFA .* sum(a,2));        % Hz, spike rate
    else
        r = relu(u_d);                            % no SFA
    end

    if n_b > 0
        p = r .* prod(b,2);                       % axonal output
    else
        p = r;                                    % no STD
    end

    %% derivatives
    if n_a > 0
        da_dt = (-a + r) ./ tau_a;
    else
        da_dt = [];
    end

    if n_b > 0
        db_dt = (1 - b) ./ tau_b - F_STD .* p ./ tau_STD;
    else
        db_dt = [];
    end

    u_d_dt = (-u_d + u_ex + M * p) ./ tau_d;

    %% load derivatives into dXdt
    dX_dt = [da_dt(:); db_dt(:); u_d_dt];

end