function [dX_dt] = SRNN(t, X, t_ex, u_ex, params)

    %% interpolate u vector
    u = interp1(t_ex, u_ex', t', 'linear')'; % u will be n x 1

    %% load parameters
    n_a = params.n_a; % number of SFA timescales per neuron
    n_b = params.n_b; % number of STD timescales per neuron

    tau_a = params.tau_a; % s, 1 x n_a, time constants of SFA
    tau_b = params.tau_b; % s, 1 x n_b, time constants of STD

    tau_d = params.tau_d; % s, scalar

    n = params.n; % n neurons
    M = params.M; % connection matrix
    
    c_SFA = params.c_SFA; % n x 1, 0 for I neurons
    F_STD = params.F_STD; % n x 1, 0 for I neurons
    tau_STD = params.tau_STD; % scalar

    %% unpack state variables using the generalized function
    % X is N_sys_eqs x 1 here.
    % unpack_SRNN_states will return:
    % a: n x n_a
    % b: n x n_b
    % u_d: n x 1
    [a, b, u_d] = unpack_SRNN_state(X, n, n_a, n_b);

    %% make dependent variables from state variables and parameters
    if n_a > 0
        r = relu(u_d - c_SFA .* sum(a,2));        % Hz, spike rate, n x 1
    else
        r = relu(u_d);                            % no SFA, n x 1
    end

    if n_b > 0
        p = r .* prod(b,2);                       % axonal output, n x 1
    else
        p = r;                                    % no STD, n x 1
    end

    %% derivatives
    % da_dt should be n x n_a
    % db_dt should be n x n_b
    % u_d_dt should be n x 1

    if n_a > 0
        % r is n x 1, tau_a is 1 x n_a. Broadcasting makes (r - a) ./ tau_a valid.
        da_dt = (r - a) ./ tau_a;
        da_dt(c_SFA==0) = 0; % enforce no adaptation for neurons where c_SFA == 0, inhibitory neurons
    else
        da_dt = [];
    end

    if n_b > 0
        % (1-b) is n x n_b. tau_b is 1 x n_b.
        % F_STD is n x 1. p is n x 1. (F_STD .* p) is n x 1.
        % Broadcasting makes (F_STD .* p) ./ tau_STD valid against (1-b)./tau_b if tau_STD is scalar.
        db_dt = (1 - b) ./ tau_b - (F_STD .* p) ./ tau_STD;
        db_dt(F_STD==0) = 0; % enforce no depression for neurons where F_STD == 0, inhibitory neurons
    else
        db_dt = [];
    end
    
    % u_d is n x 1, u is n x 1, M is n x n, p is n x 1. M*p is n x 1.
    u_d_dt = (-u_d + u + M * p) ./ tau_d;

    %% load derivatives into dXdt
    dX_dt = [da_dt(:); db_dt(:); u_d_dt];

end