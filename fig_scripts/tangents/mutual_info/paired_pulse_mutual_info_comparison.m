%--- SRNN_caller.m - Multi-condition comparison

close all
clear SRNN_NL % clear persistant variables in SRNN_NL
clearvars -except seed
clc

tic

%% 
% if exist('seed','var')
%     seed = seed+1
% else
    % seed = 105;
    seed = 103;
% end

%% Define conditions to loop over
conditions = {'none', 'sfa', 'std', 'both'};
n_conditions = length(conditions);

% Storage for results
MI_vs_delay_all = cell(n_conditions, 1);
delay_vec_all = cell(n_conditions, 1);

% Colors for plotting
colors = [0 0 1; 1 0 0; 0 0.7 0; 0.7 0 0.7]; % blue, red, green, purple

%% Loop over conditions
for i_condition = 1:n_conditions
    
    fprintf('\n=== Running condition %d/%d: %s ===\n', i_condition, n_conditions, conditions{i_condition});
    
    % Set the same seed for each condition for fair comparison
    rng(seed,'twister');
    
    % Clear persistent variables between conditions
    clear SRNN_NL
    
    %% adaptation 
    sfa_std_both_none = conditions{i_condition};
    
    tau_d = 0.1; % s, scalar
    
    %% plot saving
    save_plots = false;
    save_folder = [fullfile(fileparts(mfilename('fullpath'))) filesep 'MI_output_examples' filesep sfa_std_both_none filesep 'seed_' num2str(seed)];
    save_name = ['MI_example_' sfa_std_both_none '_tau_d_' strrep(num2str(tau_d),'.','p') '_seed_' num2str(seed)];
    
    %% plotting parameters
    sr_or_poisson = 'sr'; % sr, poisson, or sr_stacked
    include_sum_E_I_SR = false; % true or false
    
    %% Network
    n = 10; % number of neurons
    
    Lya_method = 'benettin'; % 'benettin', 'qr', 'svd', or 'none'
    use_Jacobian = false; 
    
    mean_in_out_degree = 5; % desired mean number of connections in and out
    density = mean_in_out_degree/(n-1); % each neuron can make up to n-1 connections with other neurons
    sparsity = 1-density;
    
    EI = 0.7;
    scale = 0.5/0.79782; % overall scaling factor of weights
    w.EE = scale*1.00; % E to E. Change to scale*2 for bursting
    w.EI = scale*1; % E to I connections
    w.IE = scale*1; % I to E
    w.II = scale*0.5; % I to I
    w.selfE = 0;    % self connections of E neurons
    w.selfI = 0;    % self connections of I neurons
    
    [M, EI_vec] = generate_M_no_iso(n,w,sparsity, EI);
    EI_vec = EI_vec(:); % make it a column
    [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);
    
    %% Time
    fs = 500; %Plotting sample frequency
    dt = 1/fs;
    T = [-30 600]; % [-30 2400]
    
    T_lya_1 = -15; % s, time to start Lyapunov calculation warmup
    
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
    
    % Paired-pulse stimulus parameters
    pWidth1 = 2;       % seconds, width of the first pulse
    pWidth2 = 2;       % seconds, width of the second pulse
    pAmp2 = 2.5;         % amplitude of the second pulse
    ppISI = pWidth1*1.5;         % seconds, time from pulse 1 start to pulse 2 start
    repeatISI = 15;      % seconds, time between the start of each pair
    
    % Determine the start times for each pair of pulses
    end_time = T(2);
    pair_start_times = 0:repeatISI:end_time;
    
    % Ensure the full pair, including pulse widths, fits within the simulation time
    pair_start_times = pair_start_times(pair_start_times + ppISI + pWidth2 <= end_time);
    
    % Generate random amplitudes for the first pulse of each pair
    nPairs = length(pair_start_times);
    Amp1Levels = 5; % 10
    pAmp1 = 0.5*randi([1, Amp1Levels], 1, nPairs);
    
    ch_in = 1;
    
    % Create the pulse train on channel 1
    for i = 1:nPairs
        % First pulse
        p1_start_t = pair_start_times(i);
        p1_end_t = p1_start_t + pWidth1;
        p1_indices = t >= p1_start_t & t < p1_end_t;
        u_ex(ch_in, p1_indices) = pAmp1(i);

        % Second pulse
        p2_start_t = p1_start_t + ppISI;
        p2_end_t = p2_start_t + pWidth2;
        p2_indices = t >= p2_start_t & t < p2_end_t;
        u_ex(ch_in, p2_indices) = pAmp2;
    end
    
    pulse2_start_times = pair_start_times+ppISI;
    pulse2_start_inds = max(1, min(nt, round((pulse2_start_times - T(1))*fs) + 1));
    
    DC = 0.1;
    
    % Ramp up to DC over the first 10 seconds from t=T(1)
    ramp_duration = 10; % seconds
    ramp_end_time = T(1) + ramp_duration;
    ramp_indices = t <= ramp_end_time;
    num_ramp_points = sum(ramp_indices);
    ramp_profile = linspace(0, DC, num_ramp_points);
    
    u_dc_profile = ones(1, nt) * DC;
    u_dc_profile(ramp_indices) = ramp_profile;
    u_ex = u_ex + u_dc_profile;
    
    %% parameters
    
    tau_STD = 0.5; % 0.5 scalar, time constant of synaptic depression
    
    if strcmpi(sfa_std_both_none,'both')
        n_a_E = 3; % typically 3, number of SFA timescales for E neurons
        n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
        n_b_E = 1; % typically 1 or 2, number of STD timescales for E neurons
        n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)
    elseif strcmpi(sfa_std_both_none,'sfa')
        n_a_E = 3; % typically 3, number of SFA timescales for E neurons
        n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
        n_b_E = 0; % typically 1 or 2, number of STD timescales for E neurons
        n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)
    elseif strcmpi(sfa_std_both_none,'std')
        n_a_E = 0; % typically 3, number of SFA timescales for E neurons
        n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
        n_b_E = 1; % typically 1 or 2, number of STD timescales for E neurons
        n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)
    elseif strcmpi(sfa_std_both_none,'none')
        n_a_E = 0; % typically 3, number of SFA timescales for E neurons
        n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
        n_b_E = 0; % typically 1 or 2, number of STD timescales for E neurons
        n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)
    else
        error('sfa_std_both_none must be sfa, std, both, or none')
    end
    
    % Define tau_a and tau_b for E and I neurons
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
    
    % c_SFA and F_STD remain n x 1, defining strength for *all* neurons.
    c_SFA = zeros(n, 1);
    if n_a_E > 0
        c_SFA(E_indices) = (0.5/n_a_E); % SFA strength for E neurons
    end
    if n_a_I > 0
        c_SFA(I_indices) = (0.5/n_a_I); % SFA strength for I neurons
    end
    
    F_STD = zeros(n, 1);
    if n_b_E > 0
        F_STD(E_indices) = 1; % STD for E neurons
    end
    if n_b_I > 0
        F_STD(I_indices) = 1; % STD for I neurons
    end
    
    params = package_params(n_E, n_I, E_indices, I_indices, ...
                            n_a_E, n_a_I, n_b_E, n_b_I, ...
                            tau_a_E, tau_a_I, tau_b_E, tau_b_I, ...
                            tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);
    
    params.activation_function = @(x) max(x,0);   % relu
    
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
    
    if use_Jacobian 
        ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', dt, 'InitialStep', 0.05*dt, 'Jacobian', SRNN_Jacobian_wrapper); % fast
    else
        ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e8, 'MaxStep',dt, 'InitialStep', 0.1*dt); % RelTol must be less than perturbation d0, which is 1e-3
    end
    
    SRNN_wrapper = @(tt,XX) SRNN_NL(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params to SRNN
    
    % % wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
    % solver_method = 6; % 5 is classic RK4, 6 is improved RK4
    % deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this
    % ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params))); % Pass params to get_minMaxRange
    
    % %% pick an ODE solver
    % ode_solver = @ode15s; % stiff ode solver

    % wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
    solver_method = 6; % 5 is classic RK4, 6 is improved RK4
    deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this
    ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params))); % Pass params to get_minMaxRange

    %% pick an ODE solver
    ode_solver = ode_RKn_wrapper; % fixed step RK 1, 2, or 4th order, with boundary enforcement
    % ode_solver = @ode45; % variable step
    % ode_solver = @ode4_wrapper; % basic RK4 for comparison
    % ode_solver = @ode15s; % stiff ode solver


    
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
                [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN_NL, t, u_ex, ode_solver);
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
    [a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);
    [r, p] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);
    
    % warning('adding measurement noise to r')
    r = r + (2/fs).*randn(size(r)); % measurement noise
    
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
    
    %% Make plots using the plotting function (only for first condition to avoid clutter)
    if i_condition == 1
        if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
            SRNN_tseries_figure(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, sr_or_poisson, include_sum_E_I_SR, lya_results);
        else
            SRNN_tseries_figure(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, sr_or_poisson, include_sum_E_I_SR);
        end
    end
    
    %% mutual information
    data_in = pAmp1;
    data_out = r(ch_in, pulse2_start_inds);
    
    n_bins_in = Amp1Levels;
    n_bins_out = Amp1Levels;

    %% Delayed mutual information analysis
    delay_vec = 1:10:fix(fs*pWidth2); % delays in samples from pulse 2 start
    MI_vs_delay = zeros(size(delay_vec));
    
    for i_delay = 1:length(delay_vec)
        delay_samples = delay_vec(i_delay);
        
        % Read spike rate at delayed time points during pulse 2
        readout_inds = pulse2_start_inds + delay_samples;
        
        % Make sure we don't go beyond the time series length
        valid_pairs = readout_inds <= nt;
        
        if sum(valid_pairs) > 1 % Need at least 2 data points for MI
            data_in_delayed = pAmp1(valid_pairs);
            data_out_delayed = r(ch_in, readout_inds(valid_pairs));
            
            [MI_delayed] = mutual_info_SISO(data_in_delayed, data_out_delayed, n_bins_in, n_bins_out);
            MI_vs_delay(i_delay) = MI_delayed;
        else
            MI_vs_delay(i_delay) = NaN;
        end
    end
    
    % Store results
    MI_vs_delay_all{i_condition} = MI_vs_delay;
    delay_vec_all{i_condition} = delay_vec;
    
    fprintf('Condition %s completed.\n', conditions{i_condition});
end

%% Create comparison plot
sim_dur = toc;
fprintf('\nTotal simulation duration: %.2f seconds\n', sim_dur);

% Create the comparison figure
figure('Position', [100, 100, 800, 600]);
hold on;

% Plot each condition
for i_condition = 1:n_conditions
    delay_vec = delay_vec_all{i_condition};
    MI_vs_delay = MI_vs_delay_all{i_condition};
    
    plot(delay_vec/fs, MI_vs_delay, 'LineWidth', 2.5, 'Color', colors(i_condition,:), ...
         'DisplayName', upper(conditions{i_condition}));
end

xlabel('Delay from pulse 2 start (s)')
ylabel('Mutual Information (bits)')
% title('Mutual Information vs Readout Delay - Adaptation Comparison', 'FontSize', 16)
legend('Location', 'best')
grid on
xlim([0, pWidth2])

% Add some formatting
set(gca, 'FontSize', 12);
box on;

fprintf('\n=== Final Results Summary ===\n');
for i_condition = 1:n_conditions
    fprintf('%s: Max MI = %.4f bits\n', ...
            upper(conditions{i_condition}), max(MI_vs_delay_all{i_condition}));
end
