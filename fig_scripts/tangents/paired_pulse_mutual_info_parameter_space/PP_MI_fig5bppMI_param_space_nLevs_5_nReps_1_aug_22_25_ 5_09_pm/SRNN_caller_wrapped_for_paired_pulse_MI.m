function [result] = SRNN_caller_wrapped_for_paired_pulse_MI(seed, n, EE_factor, IE_factor, EI, E_self, mean_weight, DC, mean_in_out_degree, tau_a_E_2, tau_b_E_2, tau_STD, c_SFA_factor, n_a_E, n_b_E, fs)
    % SRNN_caller_wrapped_for_paired_pulse_MI
    % Dual-stage SRNN simulation with paired-pulse stimulus and MI vs delay.
    % Matches the calling signature used by parameter space exploration.

    % Reset SRNN persistent state so stimulus changes are respected
    clear SRNN_NL
    tic

    %% Parameter validation
    n_E_expected = round(EI * n);
    if n_E_expected == 0 && c_SFA_factor > 0
        warning('EI=%.3f gives no E neurons but c_SFA_factor=%.3f>0. Setting c_SFA_factor=0', EI, c_SFA_factor);
        c_SFA_factor = 0;
    end
    if n_E_expected < 1
        error('EI * n = %.1f < 1 means no excitatory neurons. This will fail.', EI * n);
    end

    rng(seed,'twister');

    %% Network
    Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'
    use_Jacobian = false; %#ok<NASGU>

    scale = mean_weight/0.79782; % overall scaling factor of weights
    w.EE = scale*EE_factor;
    w.EI = scale*1;
    w.IE = scale*IE_factor;
    w.II = scale*.5;
    w.selfE = E_self;
    w.selfI = 0;

    density = mean_in_out_degree/(n-1);
    sparsity = 1-density;

    [M, EI_vec] = generate_M_no_iso(n, w, sparsity, EI);
    EI_vec = EI_vec(:);
    [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec); %#ok<ASGLU>

    %% Time
    dt = 1/fs;
    T = [-20 300]; % 600
    T_lya_1 = -10;

    nt = round((T(2)-T(1))*fs)+1;
    t = linspace(T(1), T(2), nt)';

    %% Paired-pulse stimulus (channel 1)
    u_ex = zeros(n, nt);

    % DC ramp
    ramp_duration = 5; % s
    ramp_end_time = T(1) + ramp_duration;
    ramp_indices = t <= ramp_end_time;
    ramp_profile = linspace(0, DC, sum(ramp_indices));
    u_dc_profile = ones(1, nt) * DC;
    u_dc_profile(ramp_indices) = ramp_profile;
    u_ex = u_ex + u_dc_profile;

    % Paired-pulse parameters (kept constant across runs for comparability)
    pWidth1 = 2;       % s
    pWidth2 = 2;       % s
    pAmp2   = 2.5;     % amplitude of pulse 2
    ppISI   = pWidth1*2.25; % *1.5 s, time from pulse1 start to pulse2 start
    repeatISI = 12;    % s, time between pair starts
    ch_in = 1;

    % Determine pair start times within [0, T(2)]
    pair_start_times = 0:repeatISI:max(0, T(2));
    % Ensure full pair fits before T(2)
    pair_start_times = pair_start_times(pair_start_times + ppISI + pWidth2 <= T(2));
    nPairs = length(pair_start_times);

    % Random first-pulse amplitudes
    Amp1Levels = 5;
    if nPairs > 0
        pAmp1 = 0.5*randi([1, Amp1Levels], 1, nPairs);
    else
        pAmp1 = [];
    end

    % Draw pulses onto u_ex for ch_in
    for iPair = 1:nPairs
        p1_start_t = pair_start_times(iPair);
        p1_end_t   = p1_start_t + pWidth1;
        p2_start_t = p1_start_t + ppISI;
        p2_end_t   = p2_start_t + pWidth2;

        p1_idx = t >= p1_start_t & t < p1_end_t;
        p2_idx = t >= p2_start_t & t < p2_end_t;
        u_ex(ch_in, p1_idx) = pAmp1(iPair);
        u_ex(ch_in, p2_idx) = pAmp2;
    end

    pulse2_start_times = pair_start_times + ppISI;
    pulse2_start_inds = max(1, min(nt, round((pulse2_start_times - T(1))*fs) + 1));

    %% Parameters
    n_a_I = 0; n_b_I = 0;
    if n_a_E > 0, tau_a_E = logspace(log10(0.3), log10(tau_a_E_2), n_a_E); else, tau_a_E = []; end
    if n_a_I > 0, tau_a_I = logspace(log10(0.3), log10(6),        n_a_I); else, tau_a_I = []; end
    if n_b_E > 0, tau_b_E = logspace(log10(0.6), log10(tau_b_E_2), n_b_E); else, tau_b_E = []; end
    if n_b_I > 0, tau_b_I = logspace(log10(0.6), log10(9),         n_b_I); else, tau_b_I = []; end

    tau_d = 0.025;

    % SFA strength for E/I groups
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

    params.activation_function = @(x) min(200, max(0, x)); % relu, satu

    %% Initial conditions
    a0_E = []; if params.n_E > 0 && params.n_a_E > 0, a0_E = zeros(params.n_E * params.n_a_E, 1); end
    a0_I = []; if params.n_I > 0 && params.n_a_I > 0, a0_I = zeros(params.n_I * params.n_a_I, 1); end
    b0_E = []; if params.n_E > 0 && params.n_b_E > 0, b0_E = ones(params.n_E * params.n_b_E, 1); end
    b0_I = []; if params.n_I > 0 && params.n_b_I > 0, b0_I = ones(params.n_I * params.n_b_I, 1); end
    u_d0 = zeros(n, 1);
    X_0 = [a0_E; a0_I; b0_E; b0_I; u_d0]; %#ok<NASGU>

    %% ODE solver
    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-7, 'MaxStep', dt, 'InitialStep', 0.01*dt);
    SRNN_wrapper = @(tt,XX) SRNN_NL(tt,XX,t,u_ex,params);

    solver_method = 6; % improved RK4
    deci = 1;          % do not decimate for Benettin
    ode_RKn_wrapper = @(odefun, tspan, y0, options) deal(tspan(:), ode_RKn_deci_bounded(odefun, tspan, y0, solver_method, false, deci, get_minMaxRange(params)));
    ode_solver = ode_RKn_wrapper;

    %% Phase 1: short stability pre-check
    LLE_phase1 = NaN; lya_results_phase1 = struct(); proceed_to_phase2 = false;
    T_phase1 = [0 5];
    idx_p1 = find(t >= T_phase1(1) & t <= T_phase1(2));
    t_p1 = t(idx_p1); u_ex_p1 = u_ex(:, idx_p1);
    [~, X_p1] = ode_solver(@(tt,XX) SRNN_NL(tt,XX,t_p1,u_ex_p1,params), t_p1, X_0, ode_options);
    T_lya_1_p1 = 0; %#ok<NASGU>
    lya_start_idx_p1 = find(t_p1 >= 0, 1, 'first');
    X_for_lya_p1 = X_p1(lya_start_idx_p1:end, :);
    t_for_lya_p1 = t_p1(lya_start_idx_p1:end);
    lya_dt_p1 = 0.5*tau_d; d0_p1 = 1e-3;
    [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya_p1, t_for_lya_p1, dt, fs, d0_p1, T_phase1, lya_dt_p1, params, ode_options, @SRNN_NL, t, u_ex, ode_solver);
    LLE_phase1 = LLE; %#ok<NASGU>
    lya_results_phase1.LLE = LLE; lya_results_phase1.local_lya = local_lya; lya_results_phase1.finite_lya = finite_lya; lya_results_phase1.t_lya = t_lya;

    [a_E_p1, a_I_p1, b_E_p1, b_I_p1, u_d_p1] = unpack_SRNN_state(X_for_lya_p1, params);
    [r_p1, ~] = compute_dependent_variables(a_E_p1, a_I_p1, b_E_p1, b_I_p1, u_d_p1, params);
    max_r_p1 = 0; 
    if ~isempty(r_p1)
        max_r_p1 = max(r_p1(:)); 
    end
    lya_results_phase1.mean_rate = mean(r_p1(:), 'omitnan');

    r_threshold = 195; % Hz, screen runaway
    proceed_to_phase2 = ~(isnan(max_r_p1) || max_r_p1 >= r_threshold);

    %% Phase 2: full run and MI if stable
    if proceed_to_phase2
        [t_ode, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options); %#ok<ASGLU>

        lya_dt = 0.5*tau_d; d0 = 1e-3;
        lya_calc_start_idx = find(t >= T_lya_1, 1, 'first');
        X_for_lya = X(lya_calc_start_idx:end, :);
        t_for_lya = t(lya_calc_start_idx:end);
        [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN_NL, t, u_ex, ode_solver);

        [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X, params);
        [r_ts, ~] = compute_dependent_variables(a_E, a_I, b_E, b_I, u_d, params);
        if ~isempty(r_ts)
            mean_rate_phase2 = mean(r_ts(:));
        else
            mean_rate_phase2 = NaN;
        end

        r_ts = r_ts + (7/fs).*randn(size(r_ts)); % some measurement noise

        lya_results = struct();
        if isfinite(LLE)
            lya_results.LLE = LLE; lya_results.local_lya = local_lya; lya_results.finite_lya = finite_lya; lya_results.t_lya = t_lya; %#ok<NASGU>
            lya_results.mean_rate = mean_rate_phase2;
        else
            lya_results = lya_results_phase1; %#ok<NASGU>
        end

        % MI vs delay (paired pulses)
        if isfinite(LLE)
            n_bins_in = Amp1Levels;
            n_bins_out = Amp1Levels;
            if isempty(pulse2_start_inds)
                MI_vs_delay_corrected = NaN; MI_vs_delay_uncorrected = NaN; delay_vec = NaN;
            else
                delay_vec = 1:10:fix(fs*pWidth2); % samples
                MI_vs_delay_corrected = zeros(size(delay_vec));
                MI_vs_delay_uncorrected = zeros(size(delay_vec));

                for iDelay = 1:length(delay_vec)
                    delay_samples = delay_vec(iDelay);
                    readout_inds = pulse2_start_inds + delay_samples;
                    valid_pairs = readout_inds <= nt;
                    if sum(valid_pairs) > 1
                        data_in_delayed = pAmp1(valid_pairs);
                        data_out_delayed = r_ts(ch_in, readout_inds(valid_pairs));
                        [MI_corrected, MI_uncorrected] = mutual_info_SISO(data_in_delayed, data_out_delayed, n_bins_in, n_bins_out);
                        MI_vs_delay_corrected(iDelay) = MI_corrected;
                        MI_vs_delay_uncorrected(iDelay) = MI_uncorrected;
                    else
                        MI_vs_delay_corrected(iDelay) = NaN;
                        MI_vs_delay_uncorrected(iDelay) = NaN;
                    end
                end
            end
            lya_results.MI_vs_delay_corrected = MI_vs_delay_corrected;
            lya_results.MI_vs_delay_uncorrected = MI_vs_delay_uncorrected;
            lya_results.delay_vec = delay_vec;
        else
            lya_results.MI_vs_delay_corrected = NaN;
            lya_results.MI_vs_delay_uncorrected = NaN;
            lya_results.delay_vec = NaN;
        end
    else
        lya_results = lya_results_phase1; %#ok<NASGU>
    end

    %% Package results
    result = struct();
    if exist('lya_results','var') && isfield(lya_results, 'LLE')
        result.LLE = lya_results.LLE;
    end
    if exist('lya_results','var') && isfield(lya_results, 'mean_rate')
        result.mean_rate = lya_results.mean_rate;
    end
    if exist('lya_results','var') && isfield(lya_results, 'MI_vs_delay_corrected')
        result.MI_vs_delay_corrected = lya_results.MI_vs_delay_corrected;
    end
    if exist('lya_results','var') && isfield(lya_results, 'MI_vs_delay_uncorrected')
        result.MI_vs_delay_uncorrected = lya_results.MI_vs_delay_uncorrected;
    end
    if exist('lya_results','var') && isfield(lya_results, 'delay_vec')
        result.delay_vec = lya_results.delay_vec;
    end

    sim_dur = toc;
    result.sim_dur = sim_dur;
    if T(2) > T(1)
        result.sim_t_dived_by_rt = sim_dur./(T(2)-T(1));
    else
        result.sim_t_dived_by_rt = NaN;
    end
end


