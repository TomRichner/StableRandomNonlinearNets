% --- Main script to estimate the Lyapunov spectrum of the Lorenz system ---
% Using QR decomposition method
clear; clc; close all;
rng(42,'twister') % For reproducibility

% --- 1. Lorenz System Parameters ---
sigma = 10.0;
rho   = 28.0;
beta  = 8/3;
params = [sigma, rho, beta]; % Store parameters in a vector
N_states = 3; % Number of states in the system

% --- 2. Algorithm Parameters ---
% Initial conditions for the fiducial trajectory
X_0 = [-4.7795941536989090537; -8.8325236589201416848; 10.425932436960884786]; % X_0 as a column vector

% --- 3. Time and sampling parameters ---
fs = 100; % Plotting sample frequency (Hz)
dt = 1/fs; % Plotting sample time (s)
T = [-100, 200]; % Plotting time interval [start, end] (s)
% T(1) to 0 is considered warmup/settling time.

% Validate time interval
if not( T(1)<=0 && 0<T(2) )
    error('T(1) must be 0 or negative, and T(2) must be positive for the LLE calculation logic.')
end

nt = round((T(2)-T(1))*fs)+1; % Number of plotting samples
t = linspace(T(1), T(2), nt)'; % Plotting time vector

lya_dt = 0.1; % Rescaling time interval for Lyapunov calculation (tau_lya) (s)

% --- 4. ODE solver options and fiducial trajectory integration ---
ode_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13, 'MaxStep', 0.1*dt, 'InitialStep', 0.05*dt);

% Integrate the fiducial trajectory
[t_ode, X_fid_traj] = ode45(@(tt,XX) lorenz_eqs(tt,XX,params), t, X_0, ode_options);

assert(all(abs(t_ode - t) < 1e-12), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
clear t_ode % t_ode is same as t

% --- 5. Lyapunov Spectrum Calculation (QR Method) -----------------
[LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = ...
    lyapunov_spectrum_qr(X_fid_traj, t, lya_dt, params, ode_options, @lorenz_jacobian_eqs, T, N_states, fs);

% Display the estimated Lyapunov Spectrum
fprintf('----------------------------------------------------\n');
fprintf('Estimated Lyapunov Spectrum (Global):\n');
for i = 1:N_states
    fprintf('  LE(%d): %f\n', i, LE_spectrum(i));
end
fprintf('Sum of exponents: %f (should be < 0 for dissipative systems like Lorenz)\n', sum(LE_spectrum));
fprintf('Kaplan-Yorke Dimension: %f\n', calculate_kaplan_yorke_dimension(LE_spectrum));
fprintf('----------------------------------------------------\n');

% --- Plot State Variables and Lyapunov Exponents ---
figure('Position', [100, 100, 800, 1000]); % Adjusted figure size
ax = gobjects(N_states + 1, 1); % Axes for states + spectrum

% Plot state variables
state_names = {'x(t)', 'y(t)', 'z(t)'};
for i = 1:N_states
    ax(i) = subplot(N_states + 1, 1, i);
    plot(t, X_fid_traj(:,i));
    xlabel('Time (t)');
    ylabel(state_names{i});
    title(['State Variable ', state_names{i}(1)]);
    grid on;
end

% Subplot for Lyapunov Spectrum
ax(N_states + 1) = subplot(N_states + 1, 1, N_states + 1);
colors = lines(N_states); % Get distinct colors for each exponent

% Plot local Lyapunov exponents
for i = 1:N_states
    plot(t_lya, local_LE_spectrum_t(:,i), '--', 'Color', colors(i,:), 'DisplayName', sprintf('Local LE(%d)', i));
    hold on;
end

% Plot finite-time Lyapunov exponents
for i = 1:N_states
    plot(t_lya, finite_LE_spectrum_t(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Finite LE(%d)', i));
end

% Plot global Lyapunov exponents as horizontal lines
for i = 1:N_states
    plot([t_lya(1) t_lya(end)], [LE_spectrum(i) LE_spectrum(i)], ':', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('Global LE(%d): %.4f', i, LE_spectrum(i)));
end
hold off;
xlabel('Time (t)');
ylabel('Lyapunov Exponent Value');
title('Lyapunov Exponent Spectrum Estimates (Local, Finite-time, Global)');
legend('show', 'Location', 'bestoutside', 'NumColumns', ceil(N_states/2));
grid on;
ylim([-max(abs(ylim)) max(abs(ylim))]); % Symmetrize y-axis if appropriate or set manually

linkaxes(ax, 'x'); % Link x-axes of all subplots

% --- 6. Plot convergence of Finite-Time Lyapunov Exponents ---
figure('Position', [950, 100, 800, 600]);
hold on;
for i = 1:N_states
    % Plot only for t >= 0 where finite-time exponents are meaningfully accumulated
    valid_indices_finite = t_lya >= 0 & ~isnan(finite_LE_spectrum_t(:,i)) & ~isinf(finite_LE_spectrum_t(:,i));
    plot(t_lya(valid_indices_finite), finite_LE_spectrum_t(valid_indices_finite,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Finite LE(%d)', i));
    % Plot global LE line for comparison
    plot([0, T(2)], [LE_spectrum(i), LE_spectrum(i)], ':', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('Global LE(%d): %.4f', i, LE_spectrum(i)));
end
hold off;
xlabel('Time (t)');
ylabel('Lyapunov Exponent Estimate');
title('Convergence of Finite-Time Lyapunov Exponents (for t \geq 0)');
legend('show', 'Location', 'best');
grid on;
ylim_current = ylim;
ylim_max_abs = max(abs(ylim_current));
ylim([-ylim_max_abs, ylim_max_abs]); % Adjust y-limits for better visualization


% --- Lorenz System Equations Function ---
function dXdt = lorenz_eqs(~, X, params)
    % Lorenz system differential equations
    sigma = params(1); rho = params(2); beta = params(3);
    x = X(1); y = X(2); z = X(3);
    dXdt = zeros(3,1);
    dXdt(1) = sigma * (y - x);
    dXdt(2) = x * (rho - z) - y;
    dXdt(3) = x * y - beta * z;
end

% --- Lorenz System Jacobian Function ---
function J = lorenz_jacobian_eqs(~, X, params)
    % Jacobian of the Lorenz system
    % Inputs:
    %   ~      : Time (not used)
    %   X      : State vector [x; y; z]
    %   params : Parameters [sigma; rho; beta]
    % Output:
    %   J      : Jacobian matrix 3x3
    sigma = params(1); rho = params(2); beta = params(3);
    x = X(1); y = X(2); z = X(3);
    J = zeros(3,3);
    J(1,1) = -sigma; J(1,2) = sigma;  J(1,3) = 0;
    J(2,1) = rho - z;J(2,2) = -1;     J(2,3) = -x;
    J(3,1) = y;      J(3,2) = x;      J(3,3) = -beta;
end

% --- Lyapunov Spectrum QR Algorithm ---
function [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya_vec] = ...
    lyapunov_spectrum_qr(X_fid_traj, t_fid_traj, lya_dt_interval, params, ode_options_main, jacobian_func_handle, T_full_interval, N_states_sys, fs_fid)
    % Calculates the Lyapunov spectrum using the QR decomposition method.
    % Integrates variational equations along a pre-computed fiducial trajectory.
    %
    % Inputs:
    %   X_fid_traj        : Fiducial trajectory (N_samples x N_states_sys)
    %   t_fid_traj        : Time vector for fiducial trajectory
    %   lya_dt_interval   : Rescaling time interval (tau_lya)
    %   params            : System parameters
    %   ode_options_main  : ODE solver options for variational equations
    %   jacobian_func_handle : Handle to the system's Jacobian function
    %   T_full_interval   : Original time interval [T_start, T_end] for global LE calculation
    %   N_states_sys      : Dimension of the system
    %   fs_fid            : Sampling frequency of the fiducial trajectory
    %
    % Outputs:
    %   LE_spectrum           : Vector of global Lyapunov exponents (N_states_sys x 1)
    %   local_LE_spectrum_t   : Matrix of local LEs over time (N_timesteps_lya x N_states_sys)
    %   finite_LE_spectrum_t  : Matrix of finite-time LEs over time (N_timesteps_lya x N_states_sys)
    %   t_lya_vec             : Time vector for Lyapunov exponent estimates

    % Create an interpolant for the fiducial trajectory for use in ODE solver
    fiducial_interpolants = cell(N_states_sys, 1);
    for i = 1:N_states_sys
        fiducial_interpolants{i} = griddedInterpolant(t_fid_traj, X_fid_traj(:,i), 'pchip');
    end

    dt_fid = 1/fs_fid;
    deci_lya = round(lya_dt_interval / dt_fid);
    if deci_lya == 0
        error('lya_dt_interval is too small compared to the fiducial trajectory sampling time, leading to zero samples per interval.');
    end
    tau_lya = dt_fid * deci_lya; 

    t_lya_indices = 1:deci_lya:length(t_fid_traj);
    t_lya_vec = t_fid_traj(t_lya_indices);

    if ~isempty(t_lya_vec) && (t_lya_vec(end) + tau_lya > t_fid_traj(end) + eps(t_fid_traj(end)))
        t_lya_vec(end) = [];
        t_lya_indices(end) = []; 
    end
    
    nt_lya = numel(t_lya_vec);
    if nt_lya == 0
        warning('No Lyapunov intervals could be formed. Check T, lya_dt_interval, and fs_fid.');
        LE_spectrum = nan(N_states_sys,1); local_LE_spectrum_t = []; finite_LE_spectrum_t = [];
        return;
    end

    Q_current = eye(N_states_sys); 
    sum_log_R_diag = zeros(N_states_sys, 1); 
    
    local_LE_spectrum_t  = zeros(nt_lya, N_states_sys);
    finite_LE_spectrum_t = zeros(nt_lya, N_states_sys);
    
    total_positive_time_accumulated = 0; 
    ode_options_var = ode_options_main;

    for k = 1:nt_lya
        t_start_segment = t_lya_vec(k);
        t_end_segment = t_start_segment + tau_lya;
        t_end_segment = min(t_end_segment, t_fid_traj(end));
        
        current_segment_duration = t_end_segment - t_start_segment;
        if current_segment_duration <= eps 
            if k > 1
                local_LE_spectrum_t(k,:) = local_LE_spectrum_t(k-1,:);
                finite_LE_spectrum_t(k,:) = finite_LE_spectrum_t(k-1,:);
            else
                local_LE_spectrum_t(k,:) = NaN;
                finite_LE_spectrum_t(k,:) = NaN;
            end
            continue;
        end
        
        t_span_ode = [t_start_segment, t_end_segment];
        Psi0_vec = reshape(Q_current, [], 1);
        
        [~, Psi_solution_vec] = ode45(@variational_eqs_ode, t_span_ode, Psi0_vec, ode_options_var);
        
        Psi_evolved_matrix = reshape(Psi_solution_vec(end,:)', [N_states_sys, N_states_sys]);
        [Q_new, R_segment] = qr(Psi_evolved_matrix);
        diag_R = diag(R_segment);
        
        log_abs_diag_R = log(abs(diag_R)); 
        valid_diag_R = abs(diag_R) > eps; 
        
        current_local_LEs = zeros(N_states_sys,1);
        current_local_LEs(valid_diag_R) = log_abs_diag_R(valid_diag_R) / current_segment_duration;
        current_local_LEs(~valid_diag_R) = -Inf; 
        local_LE_spectrum_t(k,:) = current_local_LEs';
        
        if t_start_segment >= -eps(0) 
            sum_log_R_diag(valid_diag_R) = sum_log_R_diag(valid_diag_R) + log_abs_diag_R(valid_diag_R);
            total_positive_time_accumulated = total_positive_time_accumulated + current_segment_duration;
        end
        
        if total_positive_time_accumulated > eps
            finite_LE_spectrum_t(k,:) = (sum_log_R_diag / total_positive_time_accumulated)';
        elseif k > 1 && t_start_segment >= -eps(0)
             finite_LE_spectrum_t(k,:) = finite_LE_spectrum_t(k-1,:); 
        else
            finite_LE_spectrum_t(k,:) = NaN; 
        end

        Q_current = Q_new;
    end

    if total_positive_time_accumulated > eps
        LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
    else
        warning('QR:NoPositiveTime', 'No accumulation over positive time for global LEs.');
        LE_spectrum = nan(N_states_sys,1);
    end
    
    % Nested function definition moved to the end of the parent function
    function dPsi_vec_dt = variational_eqs_ode(tt, current_Psi_vec)
        % This nested function defines the variational ODE system:
        % d(Psi)/dt = J(X_fid(t)) * Psi
        % It has access to variables from the parent function's workspace,
        % such as fiducial_interpolants, N_states_sys, jacobian_func_handle, and params.

        % Interpolate fiducial state X_fid at current time tt
        X_fid_at_tt = zeros(N_states_sys, 1);
        for state_idx_loop = 1:N_states_sys % Renamed loop variable to avoid conflict if N_states was used
            X_fid_at_tt(state_idx_loop) = fiducial_interpolants{state_idx_loop}(tt);
        end
        
        % Calculate Jacobian at X_fid(tt) using the provided function handle
        J_matrix = jacobian_func_handle(tt, X_fid_at_tt, params);
        
        % Reshape Psi_vec (input from ODE solver) to matrix form
        Psi_matrix = reshape(current_Psi_vec, [N_states_sys, N_states_sys]);
        
        % Calculate d(Psi_matrix)/dt = J * Psi
        dPsi_matrix_dt = J_matrix * Psi_matrix;
        
        % Reshape back to vector for ODE solver output
        dPsi_vec_dt = reshape(dPsi_matrix_dt, [], 1);
    end % End of nested function variational_eqs_ode

end % End of lyapunov_spectrum_qr function

% --- Function to calculate Kaplan-Yorke Dimension ---
function KY_dim = calculate_kaplan_yorke_dimension(LE_spectrum_sorted)
    % Calculates the Kaplan-Yorke dimension from a sorted Lyapunov spectrum.
    % Assumes LE_spectrum_sorted is already sorted in descending order.
    LE_spectrum_s = sort(LE_spectrum_sorted, 'descend'); % Ensure sorted
    
    KY_dim = 0;
    sum_LE = 0;
    j = 0;
    for i = 1:length(LE_spectrum_s)
        if sum_LE + LE_spectrum_s(i) >= 0
            sum_LE = sum_LE + LE_spectrum_s(i);
            j = j + 1;
        else
            % Last exponent made sum negative
            if abs(LE_spectrum_s(i)) < eps % Avoid division by zero if LE is tiny
                 KY_dim = j; % If next LE is zero, dimension is integer j
                 return;
            end
            KY_dim = j + sum_LE / abs(LE_spectrum_s(i));
            return;
        end
    end
    % If sum of all exponents is positive (e.g. expanding system or only positive LEs considered)
    KY_dim = j; 
end

