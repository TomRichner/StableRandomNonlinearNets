% plot_param_space_spot_check.m
% A script to spot-check network simulations from network_parameter_space_analysis_paired_pulse_MI_vPaper.m.
% This script loads saved parameter configurations, re-runs the simulation
% to generate full time-series data, and then creates detailed plots for visual inspection.

close all;
clear all;
clc;

%% Configuration
n_networks_rerun = 1; % Number of networks to re-run and plot for each condition
sr_or_poisson = 'sr_stacked'; % 'sr', 'poisson', or 'sr_stacked'
include_sum_E_I_SR = true;

%% Select results folder
results_base_dir = uigetdir(pwd, 'Select the parameter space results folder');
if results_base_dir == 0
    disp('No folder selected. Aborting.');
    return;
end

%% Load summary file
summary_file = fullfile(results_base_dir, 'parameter_space_analysis_summary.mat');
if ~exist(summary_file, 'file')
    error('Could not find parameter_space_analysis_summary.mat in the selected folder.');
end
summary_data = load(summary_file);
conditions = summary_data.conditions;

%% Loop through conditions
for i_cond = 1:length(conditions)
    condition_name = conditions{i_cond}.name;
    fprintf('--- Processing condition: %s ---\n', condition_name);
    
    condition_results_file = summary_data.all_conditions_summary.(condition_name).save_filename;
    if ~exist(condition_results_file, 'file')
        warning('Results file for condition %s not found. Skipping.', condition_name);
        continue;
    end
    
    fprintf('Loading results from %s...\n', condition_results_file);
    loaded_data = load(condition_results_file, 'results', 'metadata');
    results = loaded_data.results;
    metadata = loaded_data.metadata;
    num_sims = length(results);
    
    % Find successful simulations
    success_flags = cellfun(@(res) isstruct(res) && isfield(res, 'success') && res.success, results);
    successful_indices = find(success_flags);
    
    if isempty(successful_indices)
        fprintf('No successful simulations found for condition %s. Skipping.\n', condition_name);
        continue;
    end
    
    % Randomly select simulations to re-run
    num_to_run = min(n_networks_rerun, length(successful_indices));
    rand_indices = successful_indices(randperm(length(successful_indices), num_to_run));
    
    fprintf('Re-running %d simulations for condition: %s\n', num_to_run, condition_name);
    
    for i_sim = 1:num_to_run
        sim_idx = rand_indices(i_sim);
        sim_result = results{sim_idx};
        
        p = sim_result.parameters;
        seed = sim_result.seed;
        
        fprintf('  Running sim #%d (original index: %d) with seed %d...\n', i_sim, sim_idx, seed);

        % Re-run simulation to get full time series
        try
            % Debug: Print parameters being used
            fprintf('    Parameters: EE_factor=%.2f, mean_weight=%.2f, n=%d, seed=%d\n', ...
                p.EE_factor, p.mean_weight, p.n, seed);
            
            [t, X, u_ex, params_out, lya_results] = rerun_SRNN_for_tseries(p, seed);

            % Use the original time limits from simulation, but maybe show a reasonable window
            % Since T = [-30 900] in the simulation, let's show [0 50] for easier viewing
            T_plot_limits = [0 50]; % Show first 50 seconds for easier viewing, or use T for full view

            % Unpack state variables and compute dependent variables
            [a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params_out);
            [r_ts, ~] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params_out);
            
            % Get current figure count before calling SRNN_tseries_figure
            fig_count_before = length(findall(0, 'type', 'figure'));

            % Let SRNN_tseries_figure create its own figures
            SRNN_tseries_figure(t, u_ex, r_ts, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params_out, T_plot_limits, 'benettin', sr_or_poisson, include_sum_E_I_SR, lya_results);

            % Get the time series figure (the first new figure created)
            timeseries_fig_num = fig_count_before + 1;
            fig_handle = figure(timeseries_fig_num);
            
            % Build parameter string for title
            param_str_parts = {};
            for k = 1:length(metadata.params_for_grid)
                param_name = metadata.params_for_grid{k};
                if isfield(sim_result.config, param_name)
                    param_val = sim_result.config.(param_name);
                    param_str_parts{end+1} = sprintf('%s = %.2f', param_name, param_val);
                else
                    param_str_parts{end+1} = sprintf('%s = (missing)', param_name);
                end
            end
            
            % Add title to the figure created by SRNN_tseries_figure
            param_str = strjoin(param_str_parts, ', ');
            title_str = sprintf('Condition: %s (Index: %d, Seed: %d)\n%s', ...
                condition_name, sim_idx, seed, param_str);
            sgtitle(fig_handle, title_str, 'Interpreter', 'none', 'FontSize', 10);
            
            % Set figure name
            set(fig_handle, 'Name', sprintf('Condition: %s, Sim Index: %d', condition_name, sim_idx));
            
            fprintf('    Successfully created plot for condition %s, sim #%d\n', condition_name, i_sim);

        catch ME
            warning('Failed to re-run simulation for index %d. Error: %s', sim_idx, ME.message);
            fprintf('    Error details: %s\n', getReport(ME, 'extended'));
        end
    end
end

fprintf('--- Spot-check script finished. ---\n');

function [t, X, u_ex, params, lya_results] = rerun_SRNN_for_tseries(p, seed)
    % This function is a modification of SRNN_caller_wrapped_for_paired_pulse_MI
    % It runs a single simulation and returns the full time series data for plotting.
    
    clear SRNN_NL
    
    rng(seed,'twister');
    
    % Extract parameters from struct p
    n = p.n;
    EE_factor = p.EE_factor;
    IE_factor = p.IE_factor;
    EI = p.EI;
    E_self = p.E_self;
    mean_weight = p.mean_weight;
    DC = p.DC;
    mean_in_out_degree = p.mean_in_out_degree;
    tau_a_E_2 = p.tau_a_E_2;
    tau_b_E_2 = p.tau_b_E_2;
    tau_STD = p.tau_STD;
    c_SFA_factor = p.c_SFA_factor;
    n_a_E = p.n_a_E;
    n_b_E = p.n_b_E;
    fs = p.fs;

    %% Network
    scale = mean_weight/0.79782;
    w.EE = scale*EE_factor; w.EI = scale*1; w.IE = scale*IE_factor; w.II = scale*.5;
    w.selfE = E_self; w.selfI = 0;
    density = mean_in_out_degree/(n-1);
    sparsity = 1-density;
    [M, EI_vec] = generate_M_no_iso(n, w, sparsity, EI);
    EI_vec = EI_vec(:);
    [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

    %% Time
    dt = 1/fs;
    T = [-30 900];
    T_lya_1 = -15;
    nt = round((T(2)-T(1))*fs)+1;
    t = linspace(T(1), T(2), nt)';

    %% Paired-pulse stimulus (copied from SRNN_caller_wrapped_for_paired_pulse_MI)
    u_ex = zeros(n, nt);
    ramp_duration = 5; ramp_end_time = T(1) + ramp_duration;
    ramp_indices = t <= ramp_end_time;
    ramp_profile = linspace(0, DC, sum(ramp_indices));
    u_dc_profile = ones(1, nt) * DC;
    u_dc_profile(ramp_indices) = ramp_profile;
    u_ex = u_ex + u_dc_profile;

    pWidth1 = 2; pWidth2 = 1; pAmp2 = 2.5; ppISI = pWidth1*2.5;
    repeatISI = 15; ch_in = 1;
    pair_start_times = 0:repeatISI:max(0, T(2));
    pair_start_times = pair_start_times(pair_start_times + ppISI + pWidth2 <= T(2));
    nPairs = length(pair_start_times);
    Amp1Levels = 8;
    if nPairs > 0, pAmp1 = 0.5*randi([1, Amp1Levels], 1, nPairs); else, pAmp1 = []; end
    for iPair = 1:nPairs
        p1_start_t = pair_start_times(iPair); p1_end_t = p1_start_t + pWidth1;
        p2_start_t = p1_start_t + ppISI; p2_end_t = p2_start_t + pWidth2;
        p1_idx = t >= p1_start_t & t < p1_end_t;
        p2_idx = t >= p2_start_t & t < p2_end_t;
        u_ex(ch_in, p1_idx) = pAmp1(iPair);
        u_ex(ch_in, p2_idx) = pAmp2;
    end
    
    %% Parameters
    n_a_I = 0; n_b_I = 0;
    if n_a_E > 0, tau_a_E = logspace(log10(0.3), log10(tau_a_E_2), n_a_E); else, tau_a_E = []; end
    if n_a_I > 0, tau_a_I = logspace(log10(0.3), log10(6), n_a_I); else, tau_a_I = []; end
    if n_b_E > 0, tau_b_E = logspace(log10(0.6), log10(tau_b_E_2), n_b_E); else, tau_b_E = []; end
    if n_b_I > 0, tau_b_I = logspace(log10(0.6), log10(9), n_b_I); else, tau_b_I = []; end
    tau_d = 0.025;
    c_SFA = zeros(n, 1);
    if n_a_E > 0, c_SFA(EI_vec == 1) = (c_SFA_factor / n_a_E); end
    F_STD = 1 * double(EI_vec == 1);
    params = package_params(n_E, n_I, E_indices, I_indices, n_a_E, n_a_I, n_b_E, n_b_I, ...
                            tau_a_E, tau_a_I, tau_b_E, tau_b_I, tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);
    params.activation_function = @(x) min(200, max(0, x));

    %% Initial conditions
    a0_E = []; if params.n_E > 0 && params.n_a_E > 0, a0_E = zeros(params.n_E * params.n_a_E, 1); end
    a0_I = []; if params.n_I > 0 && params.n_a_I > 0, a0_I = zeros(params.n_I * params.n_a_I, 1); end
    b0_E = []; if params.n_E > 0 && params.n_b_E > 0, b0_E = ones(params.n_E * params.n_b_E, 1); end
    b0_I = []; if params.n_I > 0 && params.n_b_I > 0, b0_I = ones(params.n_I * params.n_b_I, 1); end
    u_d0 = zeros(n, 1);
    X_0 = [a0_E; a0_I; b0_E; b0_I; u_d0];

    %% ODE solver
    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-7, 'MaxStep', dt, 'InitialStep', 0.01*dt);
    SRNN_wrapper = @(tt,XX) SRNN_NL(tt,XX,t,u_ex,params);
    solver_method = 6; deci = 1;
    ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params)));
    ode_solver = ode_RKn_wrapper;

    %% Run simulation
    [~, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);
    
    %% Run LLE calculation for plotting
    lya_results = struct();
    try
        lya_dt = 0.5*tau_d; d0 = 1e-3;
        lya_calc_start_idx = find(t >= T_lya_1, 1, 'first');
        X_for_lya = X(lya_calc_start_idx:end, :);
        t_for_lya = t(lya_calc_start_idx:end);
        [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN_NL, t, u_ex, ode_solver);
        if isfinite(LLE)
            lya_results.LLE = LLE; 
            lya_results.local_lya = local_lya; 
            lya_results.finite_lya = finite_lya; 
            lya_results.t_lya = t_lya;
        end
    catch lya_ME
        fprintf('LLE calculation failed during re-run: %s\n', lya_ME.message);
    end
end
