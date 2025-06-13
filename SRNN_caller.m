% example call of SRRN()

close all
clear all % must clear all due to use of persistant variables in SRNN.m
clc

tic

%% 
seed = 7;
rng(seed,'twister');

%% Network
n = 10; % number of neurons

Lya_method = 'benettin'; % 'benettin', 'qr', 'svd', or 'none'
use_Jacobian = false;

mean_in_out_degree = 5; % desired mean number of connections in and out
density = mean_in_out_degree/(n-1); % each neuron can make up to n-1 connections with other neurons
sparsity = 1-density;

EI = 0.7;
scale = 0.5/0.79782; % overall scaling factor of weights
w.EE = scale*1; % E to E. Change to scale*2 for bursting
w.EI = scale*1; % E to I connections
w.IE = scale*1; % I to E
w.II = scale*.5; % I to I
w.selfE = 0;    % self connections of E neurons
w.selfI = 0;    % self connections of I neurons

[M, EI_vec] = generate_M(n,w,sparsity, EI);
EI_vec = EI_vec(:); % make it a column
[E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

%% Time
fs = 1000; %Plotting sample frequency
dt = 1/fs;
T = [-30 15];

T_lya_1 = -10; % s, time to start Lyapunov calculation warmup
% T_lya_1 = T(1); % s, time to start Lyapunov calculation warmup

% Validate time interval
if not( T(1)<=0 && 0<T(2) )
    error('T(1) must be 0 or negative, and T(2) must be positive for the LLE calculation logic.')
end
if ~strcmpi(Lya_method, 'none') && not(T_lya_1 < 0 && T(1) <= T_lya_1)
    warning('For Lyapunov calculations, it is recommended that T_lya_1 is negative and T(1) <= T_lya_1 to allow for a warmup period.');
end

nt = round((T(2)-T(1))*fs)+1; % Number of plotting samples
t = linspace(T(1), T(2), nt)'; % Plotting time vector

%% u_ex, external input, stimulation

u_ex = zeros(n, nt);
% % sine and square wave stim
% stim_b0 = 0.5; amp = 0.5;
% dur = 3; % duration of sine
% f_sin = 1.*ones(1,fs*dur);
% % f_sin = logspace(log10(0.5),log10(3),fs*5);
% u_ex(1,-t(1)*fs+fix(fs*6)+(1:fix(fs*dur))) = stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));
% u_ex(1,-t(1)*fs+fix(fs*1)+(1:fix(fs*dur))) = stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))');
% u_ex = u_ex*1;
% u_ex = u_ex(:,1:nt);
DC = 0.1;
u_ex = u_ex+DC;

% u_ex(:,0.2*fs:0.3*fs) = u_ex(:,0.2*fs:0.3*fs) + 0.1; % a pulse to help Lyapunov exponent to find the direction.
% u_ex(:,1:fs) = u_ex(:,1:fs)+1./fs.*randn(n,fs); % noise in the first second to help the network get off the trivial saddle node from ICs
% u_ex = u_ex+0.001./fs.*randn(n,nt); % a tiny bit of noise to help the network get off the trivial saddle node from ICs

%% add a bit of sparse noise from T(1) to min(T_lya_1+1, 0)
% if strcmpi(Lya_method,'benettin')
%     noise_indices = T(1) <= t & t <= min(T_lya_1+1, 0);
%     noise_indices = max(T_lya_1-10, T(1)) <= t & t <= min(T_lya_1+0, 0);
%     u_ex(:, noise_indices) = u_ex(:, noise_indices) + (0.0001./fs .* randn(n, sum(noise_indices))) .* (rand(1, sum(noise_indices)) < 0.05);
% end

%% parameters

tau_STD = 0.5; % scalar, time constant of synaptic depression

% Define number of timescales for E and I neurons separately
n_a_E = 3; % typically 3, number of SFA timescales for E neurons
n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
n_b_E = 1; % typically 1 or 2, number of STD timescales for E neurons
n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)

% Define tau_a and tau_b for E and I neurons
% Ensure these are empty if the corresponding n_a_X or n_b_X is 0
if n_a_E > 0
    tau_a_E = logspace(log10(0.3), log10(15), n_a_E); % s, 1 x n_a_E
else
    tau_a_E = [];
end
if n_a_I > 0
    tau_a_I = logspace(log10(0.3), log10(6), n_a_I); % s, 1 x n_a_I 
else
    tau_a_I = [];
end

if n_b_E > 0
    tau_b_E = logspace(log10(0.6), log10(9), n_b_E);  % s, 1 x n_b_E
    if n_b_E == 1 % Specific condition from original code
        tau_b_E = 4*tau_STD;
    end
else
    tau_b_E = [];
end
if n_b_I > 0
    tau_b_I = logspace(log10(0.6), log10(9), n_b_I); % s, 1 x n_b_I
    if n_b_I == 1 % Retain similar logic if ever used
        tau_b_I = 4*tau_STD;
    end
else
    tau_b_I = [];
end


tau_d = 0.025; % s, scalar

% c_SFA and F_STD remain n x 1, defining strength for *all* neurons.
% SRNN.m will use n_a_I/n_b_I to determine if states a_I/b_I exist.
if n_a_E > 0
    c_SFA = (1/n_a_E) * double(EI_vec == 1); % n x 1, Example: SFA only for E neurons
else
    c_SFA = zeros(n, 1);
end
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

%% Integrate with ODE solver

% Create Jacobian wrapper to match SRNN_wrapper signature (inclusion of params)
SRNN_Jacobian_wrapper = @(tt,XX) SRNN_Jacobian(tt,XX,params);

% ode_options = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'MaxStep', 0.5*dt, 'InitialStep', min(0.001, 0.2*dt)); % accurate
% ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep',0.05*dt,'MaxStep', 0.5*dt, 'InitialStep', 0.5*dt); % fast

% % note: ode_options are largely ignored for fix-step solvers (ode4, RKn).  For fixed-step, best to use Ralston's methods and a small step as determined by fs.  1000 Hz seems to work fine.
% % note: these option settings are important for ode45 and ode15, Benettin's method and qr method.  Need about two orders of accuracy better than the perturbation d0 in Benettin's method
if use_Jacobian 
    % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep', dt,'MaxStep', dt, 'InitialStep', dt, 'Jacobian', SRNN_Jacobian_wrapper); % fast
    ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', dt, 'InitialStep', 0.05*dt, 'Jacobian', SRNN_Jacobian_wrapper); % fast
else
    % ode_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MinStep', 0.1*dt,'MaxStep', dt, 'InitialStep', 0.5*dt); % fast
    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-7, 'MaxStep',dt, 'InitialStep', 0.01*dt); % RelTol must be less than perturbation d0, which is 1e-3
end


% ode_options = odeset('RelTol', 1e-1, 'AbsTol', 1e-1, 'MaxStep', dt, 'InitialStep', dt,'MinStep',dt); % effectively single step for when using ode45 for debug

SRNN_wrapper = @(tt,XX) SRNN(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params to SRNN

% wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
solver_method = 6; % 5 is classic RK4
deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this

ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params))); % Pass params to get_minMaxRange

%% pick an ODE solver
% ode_solver = ode_RKn_wrapper; % fixed step RK 1, 2, or 4th order, with boundary enforcement
% ode_solver = @ode45; % variable step
% ode_solver = @ode4_wrapper; % basic RK4 for comparison
ode_solver = @ode15s; % stiff ode solver

if strcmpi(Lya_method,'qr') && ~isequal(ode_solver, @ode15s)
    warning('QR method typically requires ode15s for stability. Current solver may cause issues.');
end

% Use the wrapper instead of ode15s
% [t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);

% assert(all(abs(t_ode - t) < 1e-11), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
% clear t_ode % t_ode is same as t

%% Analytic LLE Calculation
fprintf('--- Analytic LLE Calculation ---\n');
[r0_analytic, LLE_analytic] = LLE_analytic_SRNN_fcn(n, n_E, n_I, M, DC, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
if ~isnan(LLE_analytic)
    fprintf('Analytic LLE = %f\n', LLE_analytic);
    fprintf('Max analytic fixed point rate = %f Hz\n', max(r0_analytic));
end
fprintf('--------------------------------\n');

%% Two-phase LLE computation: pre-check for stability then full run

LLE_phase1 = NaN;
lya_results_phase1 = struct();
proceed_to_phase2 = false;

% --- Phase 1: Pre-check for stability ---
if ~strcmpi(Lya_method, 'none')
    fprintf('--- Phase 1: Pre-check for stability (t=0 to 1s) ---\n');

    % Set up and run Phase 1 simulation from t=0 to t=1 with original ICs
    T_phase1 = [0 1];
    t_phase1 = (T_phase1(1):dt:T_phase1(2))';
    % Use the original initial conditions X_0 for this pre-check to test for immediate divergence.
    [t_ode_p1, X_phase1] = ode_solver(SRNN_wrapper, t_phase1, X_0, ode_options);
    assert(all(abs(t_ode_p1 - t_phase1) < 1e-12), 'ODE solver did not return results exactly at the requested times for phase 1.');
    clear t_ode_p1;
    
    % Compute LLE for Phase 1
    T_lya_1_phase1 = 0;
    lya_calc_start_idx_p1 = find(t_phase1 >= T_lya_1_phase1, 1, 'first');
    X_for_lya_p1 = X_phase1(lya_calc_start_idx_p1:end, :);
    t_for_lya_p1 = t_phase1(lya_calc_start_idx_p1:end);
    
    if strcmpi(Lya_method,'qr')
        lya_dt_p1 = 4*tau_d;
    else
        lya_dt_p1 = 0.5*tau_d;
    end

    switch lower(Lya_method)
        case 'benettin'
            fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm for Phase 1...\n');
            d0_p1 = 1e-3; 
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya_p1, t_for_lya_p1, dt, fs, d0_p1, T_phase1, lya_dt_p1, params, ode_options, @SRNN, t, u_ex, ode_solver);
            LLE_phase1 = LLE;
            lya_results_phase1.LLE = LLE; lya_results_phase1.local_lya = local_lya; lya_results_phase1.finite_lya = finite_lya; lya_results_phase1.t_lya = t_lya;

        case {'qr', 'svd'}
            fprintf('Computing full Lyapunov spectrum using %s method for Phase 1...\n', upper(Lya_method));
            if strcmpi(Lya_method, 'qr')
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(X_for_lya_p1, t_for_lya_p1, lya_dt_p1, params, ode_solver, ode_options, @SRNN_Jacobian, T_phase1, N_sys_eqs, fs);
            else % svd
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_svd(X_for_lya_p1, t_for_lya_p1, lya_dt_p1, params, ode_solver, ode_options, @SRNN_Jacobian, T_phase1, N_sys_eqs, fs);
            end
            if any(isfinite(LE_spectrum))
                LLE_phase1 = max(LE_spectrum(isfinite(LE_spectrum)));
            else
                LLE_phase1 = Inf;
            end
            lya_results_phase1.LE_spectrum = LE_spectrum; lya_results_phase1.local_LE_spectrum_t = local_LE_spectrum_t; lya_results_phase1.finite_LE_spectrum_t = finite_LE_spectrum_t; lya_results_phase1.t_lya = t_lya; lya_results_phase1.N_sys_eqs = N_sys_eqs;
    end
    
    LLE_threshold = 5;
    if isnan(LLE_phase1) || LLE_phase1 < LLE_threshold
        fprintf('Phase 1 LLE = %f (< %f). Proceeding to full simulation.\n', LLE_phase1, LLE_threshold);
        proceed_to_phase2 = true;
    else
        fprintf('Phase 1 LLE = %f (>= %f). Aborting full simulation.\n', LLE_phase1, LLE_threshold);
        proceed_to_phase2 = false;
    end
else
    proceed_to_phase2 = true; % No LLE check, proceed directly
end

if proceed_to_phase2
    % --- Phase 2: Full simulation ---
    fprintf('--- Phase 2: Full simulation from T(1)=%g to T(2)=%g ---\n', T(1), T(2));
    [t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);
    assert(all(abs(t_ode - t) < 1e-11), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
    clear t_ode % t_ode is same as t
    
    % This block computes LLEs for Phase 2, and decides whether to keep them or revert to Phase 1 results
    lya_results = struct();
    phase2_LLE_is_finite = false;
    if ~strcmpi(Lya_method, 'none')
        if strcmpi(Lya_method,'qr')
            lya_dt = 4*tau_d;
        else
            lya_dt = 0.5*tau_d;
        end
        lya_calc_start_idx = find(t >= T_lya_1, 1, 'first');
        if isempty(lya_calc_start_idx)
            error('Could not find T_lya_1 in time vector t. Check T and T_lya_1 values.');
        end
        X_for_lya = X(lya_calc_start_idx:end, :);
        t_for_lya = t(lya_calc_start_idx:end);
        
        switch lower(Lya_method)
            case 'svd'
                fprintf('Computing full Lyapunov spectrum using SVD method...\n');
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_svd(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
                if any(isfinite(LE_spectrum)), phase2_LLE_is_finite = true; end
            case 'qr'
                fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
                if any(isfinite(LE_spectrum)), phase2_LLE_is_finite = true; end
            case 'benettin'
                fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
                d0 = 1e-3;
                [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex, ode_solver);
                if isfinite(LLE), phase2_LLE_is_finite = true; end
        end

        if phase2_LLE_is_finite
             if strcmpi(Lya_method, 'benettin')
                lya_results.LLE = LLE; lya_results.local_lya = local_lya; lya_results.finite_lya = finite_lya; lya_results.t_lya = t_lya;
             else % qr or svd
                lya_results.LE_spectrum = LE_spectrum; lya_results.local_LE_spectrum_t = local_LE_spectrum_t; lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; lya_results.t_lya = t_lya; lya_results.N_sys_eqs = N_sys_eqs;
             end
        else
            fprintf('Phase 2 LLE calculation was non-finite. Reverting to Phase 1 results.\n');
            lya_results = lya_results_phase1;
        end
    end
else % Did not proceed to phase 2
    fprintf('Using Phase 1 results for final output.\n');
    X = X_phase1;
    t = t_phase1;
    T = T_phase1;
    lya_results = lya_results_phase1;
end


%% Convert X to named variables and compute dependent variables for plotting and comparisons to analytic method
% Unpack using the params structure which now contains n_E, n_I, n_a_E, etc.
[a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);
% This function will need to be updated or defined to handle the new state structure
[r, p] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);

%% print the LLE analytic and compare to numerical LLE
if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
    fprintf('----------------------------------------------------\n');
    if strcmpi(Lya_method, 'benettin')
        fprintf('Estimated Largest Lyapunov Exponent (LLE): %f\n', lya_results.LLE);
    else % qr or svd
        LE_sorted = sort(lya_results.LE_spectrum,'descend');
        fprintf('Estimated Lyapunov Spectrum (Global):\n');
        for i = 1:lya_results.N_sys_eqs
            fprintf('  LE(%d): %f\n', i, LE_sorted(i));
        end
        fprintf('Sum of exponents: %f (should be < 0 for dissipative systems)\n', sum(lya_results.LE_spectrum));
        fprintf('Kaplan-Yorke Dimension: %f\n', kaplan_yorke_dim(LE_sorted));
    end
    fprintf('----------------------------------------------------\n');
    
    % --- Comparison of Analytic and Numerical Results ---
    if exist('LLE_analytic', 'var') && ~isnan(LLE_analytic)
        fprintf('--- Comparison of Results ---\n');
        % Compare LLE
        numerical_LLE = NaN;
        if strcmpi(Lya_method, 'benettin')
            numerical_LLE = lya_results.LLE;
        else % qr or svd
            if isfield(lya_results, 'LE_spectrum') && ~isempty(lya_results.LE_spectrum)
                numerical_LLE = max(lya_results.LE_spectrum(isfinite(lya_results.LE_spectrum)));
            end
        end
        if ~isnan(numerical_LLE)
            fprintf('LLE Analytic: %f  |  Numerical: %f\n', LLE_analytic, numerical_LLE);
        else
            fprintf('LLE Analytic: %f  |  Numerical: (not available)\n', LLE_analytic);
        end
        
        % Compare Fixed Point
        % Find a steady-state portion of the simulation to get numerical r0
        % Let's use the last 10% of the simulation, if t>0
        t_positive_idx = find(t>0);
        if ~isempty(t_positive_idx)
            start_idx_fp = t_positive_idx(1) + round(0.9 * numel(t_positive_idx));
            r0_numerical = mean(r(:, start_idx_fp:end), 2);
            
            fprintf('Max Fixed Point Rate (r0):\n');
            fprintf('  Analytic: %f Hz | Numerical: %f Hz\n', max(r0_analytic), max(r0_numerical));
            
            % Print norm of difference
            fp_diff_norm = norm(r0_analytic - r0_numerical);
            fprintf('Norm of difference between analytic and numerical fixed points: %f\n', fp_diff_norm);
        else
            fprintf('Could not determine numerical fixed point (no simulation time > 0).\n');
        end
        fprintf('-----------------------------\n');
    end
else
     fprintf('Skipping Lyapunov calculation - trajectory only.\n');
end

%% Make plots using the plotting function

% Call the plotting function
if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
    SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, lya_results);
else
    SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method);
end

sim_dur = toc

sim_t_dived_by_rt = sim_dur./(T(2)-T(1))