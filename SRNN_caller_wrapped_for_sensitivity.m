function [result] = SRNN_caller_wrapped_for_sensitivity(seed, n, EE_factor, IE_factor, EI, E_self, mean_weight, DC, sparsity, tau_a_E_2, tau_b_E_2, tau_STD, c_SFA_factor, n_a_E, fs)
    %SRNN_CALLER_WRAPPED_FOR_SENSITIVITY Runs a simulation of the Spiking Rate Neural Network model.
    %
    % This function wraps the SRNN simulation and is intended for use in sensitivity
    % analysis, where parameters are varied across multiple runs.
    %
    % Inputs:
    %   seed             - (scalar) Random seed for reproducibility.
    %   n                - (scalar) Number of neurons in the network.
    %   EE_factor        - (scalar) Scaling factor for excitatory-to-excitatory weights.
    %   IE_factor        - (scalar) Scaling factor for inhibitory-to-excitatory weights.
    %   EI               - (scalar) Fraction of neurons that are excitatory (between 0 and 1).
    %   E_self           - (scalar) Weight of excitatory self-connections.
    %   mean_weight      - (scalar) Mean weight of connections.
    %   DC               - (scalar) DC offset applied to external input for all neurons.
    %   sparsity         - (scalar) Fraction of connections to be removed (between 0 and 1).
    %   tau_a_E_2        - (scalar) Upper limit for the SFA time constant vector (log-spaced).
    %   tau_b_E_2        - (scalar) Upper limit for the STD time constant vector (log-spaced).
    %   tau_STD          - (scalar) Time constant for synaptic depression onset.
    %   c_SFA_factor     - (scalar) Strength of spike-frequency adaptation.
    %   n_a_E            - (scalar) Number of SFA timescales for excitatory neurons.
    %   fs               - (scalar) Sampling frequency for simulation.
    %
    % Outputs:
    %   result - (struct) A structure containing the simulation results, including:
    %            .LLE (if 'benettin')
    %            .LE_spectrum (if 'qr')
    %            .KaplanYorkeDim (if 'qr')
    %            .t (time vector)
    %            .r (firing rates)
    %            ... and other simulation data.

    tic

    %% Parameter validation
    % Check for potentially problematic parameter combinations
    n_E_expected = round(EI * n);
    if n_E_expected == 0 && c_SFA_factor > 0
        warning('EI=%.3f gives no E neurons but c_SFA_factor=%.3f>0. Setting c_SFA_factor=0', EI, c_SFA_factor);
        c_SFA_factor = 0;
    end
    if n_a_E == 0 && c_SFA_factor > 0
        warning('n_a_E=0 but c_SFA_factor=%.3f>0. SFA will have no effect', c_SFA_factor);
    end
    if n_E_expected < 1
        error('EI * n = %.1f < 1 means no excitatory neurons. This will fail.', EI * n);
    end
    if sparsity > 0.9
        warning('Very high sparsity (%.3f) may create poorly connected networks', sparsity);
    end
    if tau_a_E_2 < 0.5 && n_a_E > 0
        warning('Very fast SFA time constant (%.3f) may cause numerical issues', tau_a_E_2);
    end
    
    %% 
    rng(seed,'twister');

    %% Network
    Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'
    use_Jacobian = false;

    % sparsity is now a direct input
    scale = mean_weight/0.79782; % overall scaling factor of weights
    w.EE = scale*EE_factor; % E to E. Change to scale*2 for bursting
    w.EI = scale*1; % E to I connections
    w.IE = scale*IE_factor; % I to E
    w.II = scale*.5; % I to I
    w.selfE = E_self;    % self connections of E neurons
    w.selfI = 0;    % self connections of I neurons

    [M, EI_vec] = generate_M(n,w,sparsity, EI);
    EI_vec = EI_vec(:); % make it a column
    [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

    % Check network connectivity
    M_binary = abs(M) > 1e-10;
    if sparsity > 0.7 && n < 20
        % Simple connectivity check for small, sparse networks
        row_sums = sum(M_binary, 2);
        col_sums = sum(M_binary, 1);
        if any(row_sums == 0) || any(col_sums == 0)
            warning('Network has isolated neurons (sparsity=%.3f, n=%d)', sparsity, n);
        end
    end

    %% Time
    dt = 1/fs;
    T = [-10 11];

    % Validate time interval
    if not( T(1)<=0 && 0<T(2) )
        error('T(1) must be 0 or negative, and T(2) must be positive for the LLE calculation logic.')
    end

    nt = round((T(2)-T(1))*fs)+1; % Number of plotting samples
    t = linspace(T(1), T(2), nt)'; % Plotting time vector

    %% u_ex, external input, stimulation

    u_ex = zeros(n, nt);
    % sine and square wave stim
    stim_b0 = 0.5; amp = 0.5;
    dur = 3; % duration of sine
    f_sin = 1.*ones(1,fs*dur);
    % f_sin = logspace(log10(0.5),log10(3),fs*5);
    u_ex(1,-t(1)*fs+fix(fs*6)+(1:fix(fs*dur))) = stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));
    u_ex(1,-t(1)*fs+fix(fs*1)+(1:fix(fs*dur))) = stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))');
    u_ex = u_ex*1;
    u_ex = u_ex(:,1:nt);
    u_ex = u_ex+DC;

    u_ex(:,0.2*fs:0.3*fs) = u_ex(:,0.2*fs:0.3*fs) + 0.1; % a pulse to help Lyapunov exponent to find the direction.
    u_ex(:,1:fs) = u_ex(:,1:fs)+1./fs.*randn(n,fs); % noise in the first second to help the network get off the trivial saddle node from ICs
    u_ex = u_ex+0.001./fs.*randn(n,nt); % a tiny bit of noise to help the network get off the trivial saddle node from ICs

    % noise_density = 0.02; % Define the density for sparse noise application
    % u_ex = u_ex + (0.001./fs .* randn(n, nt)) .* (rand(1, nt) < noise_density); % Apply sparse noise (density ~0.02) to help the network get off the trivial saddle node from ICs


    %% parameters

    % Define number of timescales for E and I neurons separately
    n_a_I = 0; % number of SFA timescales for I neurons (typically 0)
    n_b_E = 1; % number of STD timescales for E neurons
    n_b_I = 0; % number of STD timescales for I neurons (typically 0)

    % Define tau_a and tau_b for E and I neurons
    % Ensure these are empty if the corresponding n_a_X or n_b_X is 0

    if n_a_E > 0
        tau_a_E = logspace(log10(0.3), log10(tau_a_E_2), n_a_E); % s, 1 x n_a_E
    else
        tau_a_E = [];
    end
    if n_a_I > 0
        tau_a_I = logspace(log10(0.3), log10(6), n_a_I); % s, 1 x n_a_I 
    else
        tau_a_I = [];
    end

    if n_b_E > 0
        tau_b_E = logspace(log10(0.6), log10(tau_b_E_2), n_b_E);  % s, 1 x n_b_E
    else
        tau_b_E = [];
    end

    if n_b_I > 0
        tau_b_I = logspace(log10(0.6), log10(9), n_b_I); % s, 1 x n_b_I
    else
        tau_b_I = [];
    end


    tau_d = 0.025; % s, scalar

    % c_SFA and F_STD remain n x 1, defining strength for *all* neurons.
    % SRNN.m will use n_a_I/n_b_I to determine if states a_I/b_I exist.
    c_SFA = c_SFA_factor * double(EI_vec == 1); % n x 1, Example: SFA only for E neurons
    % c_SFA(I_indices) = 0; % Explicitly set to 0 for I if desired, or rely on n_a_I = 0
    F_STD = 1 * double(EI_vec == 1); % n x 1, Example: STD only for E neurons
    % F_STD(I_indices) = 0; % Explicitly set to 0 for I if desired, or rely on n_b_I = 0


    params = package_params(n_E, n_I, E_indices, I_indices, ...
                            n_a_E, n_a_I, n_b_E, n_b_I, ...
                            tau_a_E, tau_a_I, tau_b_E, tau_b_I, ...
                            tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);

    %% Initial Conditions
    a0_E = [];
    if params.n_E > 0 && params.n_a_E > 0
        a0_E = zeros(params.n_E * params.n_a_E, 1);
    end

    a0_I = [];
    if params.n_I > 0 && params.n_a_I > 0
        a0_I = zeros(params.n_I * params.n_a_I, 1);
    end

    b0_E = [];
    if params.n_E > 0 && params.n_b_E > 0
        b0_E = ones(params.n_E * params.n_b_E, 1);
    end

    b0_I = [];
    if params.n_I > 0 && params.n_b_I > 0
        b0_I = ones(params.n_I * params.n_b_I, 1);
    end

    u_d0 = zeros(n, 1);

    X_0 = [a0_E; a0_I; b0_E; b0_I; u_d0];

    N_sys_eqs = size(X_0,1); % Number of system equations / states

    % Memory usage check
    estimated_memory_mb = N_sys_eqs * nt * 8 / 1e6; % rough estimate in MB
    if estimated_memory_mb > 500 % 500 MB threshold
        warning('Simulation may use ~%.1f MB memory (N_sys_eqs=%d, nt=%d)', estimated_memory_mb, N_sys_eqs, nt);
    end

    %% Integrate with ODE solver

    % Create Jacobian wrapper to match SRNN_wrapper signature (inclusion of params)
    SRNN_Jacobian_wrapper = @(tt,XX) SRNN_Jacobian(tt,XX,params);

    % ode_options = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'MaxStep', 0.5*dt, 'InitialStep', min(0.001, 0.2*dt)); % accurate
    % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep',0.05*dt,'MaxStep', 0.5*dt, 'InitialStep', 0.5*dt); % fast

    % % note: ode_options are largely ignored for fix-step solvers (ode4, RKn).  For fixed-step, best to use Ralston's methods and a small step as determined by fs.  1000 Hz seems to work fine.
    % % note: these option settings are important for ode45 and ode15, Benettin's method and qr method.  Need about two orders of accuracy better than the perturbation d0 in Benettin's method
    if use_Jacobian 
        % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep', dt,'MaxStep', dt, 'InitialStep', dt, 'Jacobian', SRNN_Jacobian_wrapper); % fast
        ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', dt, 'InitialStep', 0.1*dt, 'Jacobian', SRNN_Jacobian_wrapper); % fast
    else
        % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep', 0.1*dt,'MaxStep', dt, 'InitialStep', 0.5*dt); % fast
        ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep',dt, 'InitialStep', 0.1*dt); % RelTol must be less than perturbation d0, which is 1e-3
    end


    % ode_options = odeset('RelTol', 1e-1, 'AbsTol', 1e-1, 'MaxStep', dt, 'InitialStep', dt,'MinStep',dt); % effectively single step for when using ode45 for debug

    SRNN_wrapper = @(tt,XX) SRNN(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params to SRNN

    % wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
    solver_method = 3; % 5 is classic RK4
    deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this
    ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params))); % Pass params to get_minMaxRange

    %% pick an ODE solver
    % ode_solver = ode_RKn_wrapper; % fixed step RK 1, 2, or 4th order, with boundary enforcement
    % ode_solver = @ode45; % variable step
    % ode_solver = @ode4_wrapper; % basic RK4 for comparison
    ode_solver = @ode15s; % stiff ode solver

    % Use the wrapper instead of ode15s
    [t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);

    assert(all(abs(t_ode - t) < 1e-11), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
    clear t_ode % t_ode is same as t

    %% comput LLE or Lyapunov spectrum

    if strcmpi(Lya_method,'qr')
        lya_dt = 0.05; % longer better for qr?
    else
        lya_dt = 0.005; % 0.005 is good for Benettin.  Rescaling time interval for Lyapunov calculation (tau_lya) (s)
    end

    switch lower(Lya_method)
        case 'qr'
            fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
            
            % Ensure SRNN_jacobian_eqs is defined elsewhere or this will error.
            % Using N_sys_eqs for the number of states.
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = ...
                lyapunov_spectrum_qr(X, t, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);

            LE_sorted = sort(LE_spectrum,'descend');
            KY_dim = kaplan_yorke_dim(LE_sorted);
            % Display the estimated Lyapunov Spectrum
            fprintf('----------------------------------------------------\n');
            fprintf('Estimated Lyapunov Spectrum (Global):\n');
            for i = 1:N_sys_eqs
                fprintf('  LE(%d): %f\n', i, LE_sorted(i));
            end
            fprintf('Sum of exponents: %f (should be < 0 for dissipative systems)\n', sum(LE_spectrum));
            fprintf('Kaplan-Yorke Dimension: %f\n', KY_dim);
            fprintf('----------------------------------------------------\n');
            
        case 'benettin'
            fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
            
            d0 = 1e-3; % Initial separation magnitude for Benettin's algorithm
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex, ode_solver);

            fprintf('----------------------------------------------------\n');
            fprintf('Estimated Largest Lyapunov Exponent (LLE): %f\n', LLE);
            fprintf('----------------------------------------------------\n');
            
        case 'none'
            fprintf('Skipping Lyapunov calculation - trajectory only.\n');
            
        otherwise
            error('Unknown method: %s. Choose ''qr'', ''benettin'', or ''none''.', method);
    end

    %% Convert X to named variables
    % Unpack using the params structure which now contains n_E, n_I, n_a_E, etc.
    [a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);

    % compute dependent variables r and p
    %% Compute dependent variables r and p using a subfunction
    % This function will need to be updated or defined to handle the new state structure
    [r, p] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);

    %% Make plots using the plotting function

    % Prepare Lyapunov results structure if needed
    lya_results = struct();
    if ~strcmpi(Lya_method, 'none')
        if strcmpi(Lya_method, 'benettin')
            if exist('LLE', 'var'), lya_results.LLE = LLE; end
            if exist('local_lya', 'var'), lya_results.local_lya = local_lya; end
            if exist('finite_lya', 'var'), lya_results.finite_lya = finite_lya; end
            if exist('t_lya', 'var'), lya_results.t_lya = t_lya; end
        elseif strcmpi(Lya_method, 'qr')
            if exist('LE_spectrum', 'var'), lya_results.LE_spectrum = LE_spectrum; end
            if exist('local_LE_spectrum_t', 'var'), lya_results.local_LE_spectrum_t = local_LE_spectrum_t; end
            if exist('finite_LE_spectrum_t', 'var'), lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; end
            if exist('t_lya', 'var'), lya_results.t_lya = t_lya; end
            if exist('N_sys_eqs', 'var'), lya_results.N_sys_eqs = N_sys_eqs; end
            if exist('KY_dim', 'var'), lya_results.KaplanYorkeDim = KY_dim; end
        end
    end

    % For sensitivity analysis, plotting is commented out.
    % if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
    %     SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, lya_results);
    % else
    %     SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method);
    % end

    result = lya_results;
    result.t = t;
    result.r = r;
    result.p = p;
    result.u_ex = u_ex;
    result.params = params;
    result.X = X;
    result.a_E_ts = a_E_ts;
    result.a_I_ts = a_I_ts;
    result.b_E_ts = b_E_ts;
    result.b_I_ts = b_I_ts;
    result.u_d_ts = u_d_ts;

    sim_dur = toc;

    result.sim_dur = sim_dur;
    result.sim_t_dived_by_rt = sim_dur./(T(2)-T(1));

end