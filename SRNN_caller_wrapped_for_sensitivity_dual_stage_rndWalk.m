function [result] = SRNN_caller_wrapped_for_sensitivity_dual_stage_rndWalk(seed, n, EE_factor, IE_factor, EI, E_self, mean_weight, DC, mean_in_out_degree, tau_a_E_2, tau_b_E_2, tau_STD, c_SFA_factor, n_a_E, n_b_E, fs)
    %SRNN_CALLER_WRAPPED_FOR_SENSITIVITY_DUAL_STAGE Runs a simulation of the Spiking Rate Neural Network model.
    %
    % This function wraps the SRNN simulation and is intended for use in sensitivity
    % analysis. It uses a dual-stage approach from SRNN_caller.m to first check for 
    % stability in a short simulation before running the full simulation, which can 
    % save significant time by aborting unstable runs early.

    %–– reset the SRNN interpolant so any new u_ex values get picked up
    clear SRNN % clear persistant variables in SRNN.m in case u_ex changes 
    tic

    %% Parameter validation (from SRNN_caller_wrapped_for_sensitivity.m)
    n_E_expected = round(EI * n);
    if n_E_expected == 0 && c_SFA_factor > 0
        warning('EI=%.3f gives no E neurons but c_SFA_factor=%.3f>0. Setting c_SFA_factor=0', EI, c_SFA_factor);
        c_SFA_factor = 0;
    end
    
    if n_E_expected < 1
        error('EI * n = %.1f < 1 means no excitatory neurons. This will fail.', EI * n);
    end
    
    %% 
    rng(seed,'twister');

    %% Network
    Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'
    use_Jacobian = false;

    scale = mean_weight/0.79782; % overall scaling factor of weights
    w.EE = scale*EE_factor; 
    w.EI = scale*1; 
    w.IE = scale*IE_factor;
    w.II = scale*.5; 
    w.selfE = E_self;
    w.selfI = 0;

    density = mean_in_out_degree/(n-1);
    sparsity = 1-density;

    [M, EI_vec] = generate_M_no_iso(n,w,sparsity, EI);
    EI_vec = EI_vec(:);
    [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

    M_binary = abs(M) > 1e-10;
    if sparsity > 0.7 && n < 20
        row_sums = sum(M_binary, 2);
        col_sums = sum(M_binary, 1);
        if any(row_sums == 0) || any(col_sums == 0)
            warning('Network has isolated neurons (sparsity=%.3f, n=%d)', sparsity, n);
        end
    end

    %% Time
    dt = 1/fs;
    T = [-20 30];
    T_lya_1 = -10;

    nt = round((T(2)-T(1))*fs)+1;
    t = linspace(T(1), T(2), nt)';

    %% External Input (u_ex)
    u_ex = zeros(n, nt);
    
    % Ramp up to DC over the first 3 seconds from t=T(1)
    ramp_duration = 3; % seconds
    ramp_end_time = T(1) + ramp_duration;
    ramp_indices = t <= ramp_end_time;
    num_ramp_points = sum(ramp_indices);
    ramp_profile = linspace(0, DC, num_ramp_points);

    u_dc_profile = ones(1, nt) * DC;
    u_dc_profile(ramp_indices) = ramp_profile;
    u_ex = u_ex + u_dc_profile;

    % add random walk stimulus to 3 neurons
    [bH,aH] = butter(1,0.1/(fs/2),'high'); 
    [bL,aL] = butter(3,10/(fs/2),'low');
    u_ex(1:3,:) = u_ex(1:3,:)+filter(bL,aL,filter(bH,aH,cumsum(10./fs.*randn(3,nt),2),[],2),[],2);


    %% add a bit of sparse noise
    % if strcmpi(Lya_method,'benettin')
    %     noise_indices = max(T_lya_1-20, T(1)) <= t & t <= min(T_lya_1-5, 0);
    %     u_ex(:, noise_indices) = u_ex(:, noise_indices) + (0.0001./fs .* randn(n, sum(noise_indices))) .* (rand(1, sum(noise_indices)) < 0.05);
    % end

    %% parameters
    n_a_I = 0; n_b_I = 0;

    if n_a_E > 0, tau_a_E = logspace(log10(0.3), log10(tau_a_E_2), n_a_E); else, tau_a_E = []; end
    if n_a_I > 0, tau_a_I = logspace(log10(0.3), log10(6), n_a_I); else, tau_a_I = []; end
    if n_b_E > 0, tau_b_E = logspace(log10(0.6), log10(tau_b_E_2), n_b_E); else, tau_b_E = []; end
    if n_b_I > 0, tau_b_I = logspace(log10(0.6), log10(9), n_b_I); else, tau_b_I = []; end
    
    tau_d = 0.025;
    % SFA strength: assign for E and I neurons separately, depending on n_a_E and n_a_I
    c_SFA = zeros(n, 1);
    if n_a_E > 0
        c_SFA(EI_vec == 1) = (c_SFA_factor / n_a_E);
    end
    if n_a_I > 0
        c_SFA(EI_vec == -1) = (c_SFA_factor / n_a_I);
    end
    F_STD = 1 * double(EI_vec == 1);

    params = package_params(n_E, n_I, E_indices, I_indices, n_a_E, n_a_I, n_b_E, n_b_I, ...
                        tau_a_E, tau_a_I, tau_b_E, tau_b_I, tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);

    %% Initial Conditions
    a0_E = []; if params.n_E > 0 && params.n_a_E > 0, a0_E = zeros(params.n_E * params.n_a_E, 1); end
    a0_I = []; if params.n_I > 0 && params.n_a_I > 0, a0_I = zeros(params.n_I * params.n_a_I, 1); end
    b0_E = []; if params.n_E > 0 && params.n_b_E > 0, b0_E = ones(params.n_E * params.n_b_E, 1); end
    b0_I = []; if params.n_I > 0 && params.n_b_I > 0, b0_I = ones(params.n_I * params.n_b_I, 1); end
    u_d0 = zeros(n, 1);
    X_0 = [a0_E; a0_I; b0_E; b0_I; u_d0];
    N_sys_eqs = size(X_0,1);

    %% ODE Solver Setup
    % ode_options = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'MaxStep', 0.5*dt, 'InitialStep', min(0.001, 0.2*dt)); % accurate
    % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep',0.05*dt,'MaxStep', 0.5*dt, 'InitialStep', 0.5*dt); % fast
    % % note: these option settings are important for ode45 and ode15, Benettin's method and qr method.  Need about two orders of accuracy better than the perturbation d0 in Benettin's method
    % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep', 0.1*dt,'MaxStep', dt, 'InitialStep', 0.5*dt); % fast
    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-7, 'MaxStep',dt, 'InitialStep', 0.01*dt); % RelTol must be less than perturbation d0, which is 1e-3
    % ode_options = odeset('RelTol', 1e-1, 'AbsTol', 1e-1, 'MaxStep', dt, 'InitialStep', dt,'MinStep',dt); % effectively single step for when using ode45 for debug
    
    SRNN_wrapper = @(tt,XX) SRNN(tt,XX,t,u_ex,params);
    
    % % wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
    solver_method = 6; % 5 is classic RK4
    deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this
    ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params))); % Pass params to get_minMaxRange
    ode_solver = ode_RKn_wrapper; % fixed step RK 1, 2, or 4th order, with boundary enforcement
    % ode_solver = @ode45; % variable step
    % ode_solver = @ode4_wrapper; % basic RK4 for comparison
    % ode_solver = @ode15s;

    %% Two-phase LLE computation
    LLE_phase1 = NaN;
    lya_results_phase1 = struct();
    proceed_to_phase2 = false;

    % --- Phase 1: Pre-check for stability ---
    T_phase1 = [0 5];
    % find which samples of the full time vector `t` lie in [0,1]
    idx_phase1 = find(t >= T_phase1(1) & t <= T_phase1(2));
    % extract the matching time vector and input slice
    t_phase1     = t(idx_phase1);
    u_ex_phase1  = u_ex(:, idx_phase1);
    % now integrate using only the first‐second of u_ex
    [~, X_phase1] = ode_solver(@(tt,XX) SRNN(tt,XX,t_phase1,u_ex_phase1,params), t_phase1, X_0, ode_options);
    
    T_lya_1_phase1 = 0;
    lya_calc_start_idx_p1 = find(t_phase1 >= T_lya_1_phase1, 1, 'first');
    X_for_lya_p1 = X_phase1(lya_calc_start_idx_p1:end, :);
    t_for_lya_p1 = t_phase1(lya_calc_start_idx_p1:end);
    
    lya_dt_p1 = 0.5*tau_d;
    d0_p1 = 1e-3; 
    [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya_p1, t_for_lya_p1, dt, fs, d0_p1, T_phase1, lya_dt_p1, params, ode_options, @SRNN, t, u_ex, ode_solver);
    LLE_phase1 = LLE;
    lya_results_phase1.LLE = LLE; lya_results_phase1.local_lya = local_lya; lya_results_phase1.finite_lya = finite_lya; lya_results_phase1.t_lya = t_lya;

    [a_E_p1, a_I_p1, b_E_p1, b_I_p1, u_d_p1] = unpack_SRNN_state(X_for_lya_p1, params);
    [r_p1, ~] = compute_dependent_variables(a_E_p1, a_I_p1, b_E_p1, b_I_p1, u_d_p1, params);
    
    max_r_phase1 = 0;
    if ~isempty(r_p1)
        max_r_phase1 = max(r_p1(:));
        mean_rate_phase1 = mean(r_p1(:));
    else
        mean_rate_phase1 = NaN;
    end
    lya_results_phase1.mean_rate = mean_rate_phase1;

    r_threshold = 1000; % Hz
    if isnan(max_r_phase1) || max_r_phase1 < r_threshold
        proceed_to_phase2 = true;
    else
        proceed_to_phase2 = false;
    end

    if proceed_to_phase2
        % --- Phase 2: Full simulation ---
        [t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);
        
        lya_results = struct();
        phase2_LLE_is_finite = false;
        lya_dt = 0.5*tau_d;
        lya_calc_start_idx = find(t >= T_lya_1, 1, 'first');
        X_for_lya = X(lya_calc_start_idx:end, :);
        t_for_lya = t(lya_calc_start_idx:end);
        
        d0 = 1e-3;
        [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex, ode_solver);
        if isfinite(LLE), phase2_LLE_is_finite = true; end

        [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X_for_lya, params);
        [r_ts, ~] = compute_dependent_variables(a_E, a_I, b_E, b_I, u_d, params);
        if ~isempty(r_ts)
            mean_rate_phase2 = mean(r_ts(:));
        else
            mean_rate_phase2 = NaN;
        end
        
        if phase2_LLE_is_finite
            lya_results.LLE = LLE; lya_results.local_lya = local_lya; lya_results.finite_lya = finite_lya; lya_results.t_lya = t_lya;
            lya_results.mean_rate = mean_rate_phase2;
        else
            lya_results = lya_results_phase1;
        end
    else % Did not proceed to phase 2
        lya_results = lya_results_phase1;
    end

    %% Package results for sensitivity analysis
    result = struct();
    if isfield(lya_results, 'LLE')
        result.LLE = lya_results.LLE;
    end
    if isfield(lya_results, 'mean_rate')
        result.mean_rate = lya_results.mean_rate;
    end

    sim_dur = toc;
    result.sim_dur = sim_dur;
    if T(2) > T(1)
      result.sim_t_dived_by_rt = sim_dur./(T(2)-T(1));
    else
      result.sim_t_dived_by_rt = NaN;
    end

end
