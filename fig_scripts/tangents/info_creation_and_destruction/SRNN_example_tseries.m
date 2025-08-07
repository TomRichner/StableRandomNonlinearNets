%--- SRNN_caller.m

close all
clear all % must clear all due to use of persistant variables in SRNN.m
clc

tic

%% 
seed = 42;
rng(seed,'twister');

%% plotting parameters
sr_or_poisson = 'sr_stacked'; % sr, poisson, or sr_stacked
include_sum_E_I_SR = false; % true or false

%% Network
n = 25; % number of neurons

Lya_method = 'benettin'; % 'benettin', 'qr', 'svd', or 'none'
use_Jacobian = false; 

mean_in_out_degree = 4; % desired mean number of connections in and out
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

[M, EI_vec] = generate_M_no_iso(n,w,sparsity, EI);
EI_vec = EI_vec(:); % make it a column
[E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

%% Time
fs = 1000; %Plotting sample frequency
dt = 1/fs;
T = [-30 20];

T_lya_1 = -20; % s, time to start Lyapunov calculation warmup
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
% sine and square wave stim
stim_b0 = 0.5; amp = 0.5;
dur = 3*1; % duration of sine
f_sin = 1.*ones(1,fs*dur);
% u_ex(1:3,-t(1)*fs+fix(fs*30)+(1:fix(fs*dur))) = 1*[1;2;3]*(stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))')));
% u_ex(1:3,-t(1)*fs+fix(fs*5)+(1:fix(fs*dur))) = 1*[1;2;3]*(stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));
u_ex(1:3,-t(1)*fs+fix(fs*8)+(1:fix(fs*dur))) = 0.5*[1;2;3]*(stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))')));
u_ex(1:3,-t(1)*fs+fix(fs*3)+(1:fix(fs*dur))) = 0.5*[1;2;3]*(stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));

u_ex = 1.5*u_ex;
% u_ex = u_ex(:,1:nt);
DC = 0.01;
% Ramp up to DC over the first 5 seconds from t=T(1)
ramp_duration = 10; % seconds
ramp_end_time = T(1) + ramp_duration;
ramp_indices = t <= ramp_end_time;
num_ramp_points = sum(ramp_indices);
ramp_profile = linspace(0, DC, num_ramp_points);

u_dc_profile = ones(1, nt) * DC;
u_dc_profile(ramp_indices) = ramp_profile;
u_ex = u_ex + u_dc_profile;

% u_ex = u_ex + (1./fs).*randn(size(u_ex)); % add some noise

%% parameters

tau_STD = 0.5; % scalar, time constant of synaptic depression

% Define number of timescales for E and I neurons separately
n_a_E = 1; % typically 3, number of SFA timescales for E neurons
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
    tau_a_I = logspace(log10(0.3), log10(15), n_a_I); % s, 1 x n_a_I 
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
    c_SFA = (0.5/n_a_E) * double(EI_vec == 1); % n x 1, Example: SFA only for E neurons
else
    c_SFA = zeros(n, 1);
end

F_STD = 1 * double(EI_vec == 1); % n x 1, Example: STD only for E neurons

params = package_params(n_E, n_I, E_indices, I_indices, ...
                        n_a_E, n_a_I, n_b_E, n_b_I, ...
                        tau_a_E, tau_a_I, tau_b_E, tau_b_I, ...
                        tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);

params.activation_function = @(x) max(min(x+0.01.*x.^1, 2),0);   % nonlinear activation function
% params.activation_function = @(x) max(min(x,10),0);   % nonlinear activation function
% params.activation_function = @(x) max(x,0);   % nonlinear activation function

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
    ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep',dt, 'InitialStep', 0.1*dt); % RelTol must be less than perturbation d0, which is 1e-3
end

% ode_options = odeset('RelTol', 1e-1, 'AbsTol', 1e-1, 'MaxStep', dt, 'InitialStep', dt,'MinStep',dt); % effectively single step for when using ode45 for debug

SRNN_wrapper = @(tt,XX) SRNN_NL(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params to SRNN

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

%% Run Simulation and LLE Computation
fprintf('--- Running simulation from T(1)=%g to T(2)=%g ---\n', T(1), T(2));
[t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);
assert(all(abs(t_ode - t) < 1e-12), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
clear t_ode % t_ode is same as t

lya_results = struct();
if ~strcmpi(Lya_method, 'none')
    fprintf('--- LLE Calculation (%s method) ---\n', upper(Lya_method));
    
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
        case {'qr', 'svd'}
            if strcmpi(Lya_method, 'qr')
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
            else % svd
                [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_svd(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
            end
            if ~any(isfinite(LE_spectrum))
                warning('LLE calculation resulted in non-finite values.');
            end
            lya_results.LE_spectrum = LE_spectrum; 
            lya_results.local_LE_spectrum_t = local_LE_spectrum_t; 
            lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; 
            lya_results.t_lya = t_lya; 
            lya_results.N_sys_eqs = N_sys_eqs;

        case 'benettin'
            d0 = 1e-3;
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex, ode_solver);
            if ~isfinite(LLE)
                warning('Benettin LLE calculation resulted in a non-finite value.');
            end
            lya_results.LLE = LLE; 
            lya_results.local_lya = local_lya; 
            lya_results.finite_lya = finite_lya; 
            lya_results.t_lya = t_lya;
    end
end


%% Convert X to named variables and compute dependent variables for plotting and comparisons to analytic method
% Unpack using the params structure which now contains n_E, n_I, n_a_E, etc.
[a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);
% This function will need to be updated or defined to handle the new state structure
[r, p] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);

%% print the LLE
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
else
     fprintf('Skipping Lyapunov calculation - trajectory only.\n');
end

%% Make plots using the plotting function

% Call the plotting function
% if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
%     SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, lya_results);
% else
%     SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method);
% end

if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
    SRNN_tseries_figure(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, sr_or_poisson, include_sum_E_I_SR, lya_results);
else
    SRNN_tseries_figure(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, sr_or_poisson, include_sum_E_I_SR);
end

figure(1)
% xlim([0 T(2)])
% subplot(6,1,6)
% ylim([-50 15])

% add subplots

sim_dur = toc

sim_t_dived_by_rt = sim_dur./(T(2)-T(1))


%% Add letters to plot
fig1 = figure(1);
% Get axes, excluding legends/colorbars, to determine how many letters.
axes_in_fig = findall(fig1, 'type', 'axes', '-not', {'Tag', 'legend'}, '-not', {'Tag', 'Colorbar'});
num_subplots = numel(axes_in_fig);

if num_subplots > 0
    if num_subplots <= 26
        letters = arrayfun(@(x) sprintf('(%c)', x), 'a':char('a'+num_subplots-1), 'UniformOutput', false);
        AddLetters2Plots(fig1, letters, 'FontSize', 18, 'FontWeight', 'Normal', 'HShift', -0.025, 'VShift', -0.032, 'Location', 'NorthWest');
    else
        warning('More than 26 subplots in Figure 1, not adding letters.');
    end
end


%% save the figs
warning('not saving figs')
% disp('Saving figures...');
% save_folder = fullfile(fileparts(mfilename('fullpath')), 'v4_output_figs_sr_stacked');
% save_name = ['example_sr_stacked_tseries_seed' num2str(seed)];
% save_some_figs_to_folder_2(save_folder, save_name, [], []);
% disp(['Figures saved to ' save_folder]);

% save_data_figs_mfiles(input_folder_path, output_folder_path, folder_name, note_string, save_mat_files, save_m_files, save_open_figs, varargin) % could also save current mfile for reproducibility

