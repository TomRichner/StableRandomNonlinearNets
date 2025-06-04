function params = package_params(n_a, n_b, tau_a, tau_b, tau_d, n, M, c_SFA, F_STD, tau_STD)
    % these params are sent into the ode solver via an anonymous function, nice to package them
    params.n_a = n_a;
    params.n_b = n_b;
    params.tau_a = tau_a;
    params.tau_b = tau_b;
    params.tau_d = tau_d;
    params.n = n;
    params.M = M;
    params.c_SFA = c_SFA;
    params.F_STD = F_STD;
    params.tau_STD = tau_STD;
end