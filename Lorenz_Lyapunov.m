% --- Main script to estimate the largest Lyapunov exponent of the Lorenz system ---
% Using Benettin's algorithm

clear; clc; close all;

rng(42,'twister')

% --- 1. Lorenz System Parameters ---
sigma = 10.0;
rho   = 28.0;
beta  = 8/3;
params = [sigma, rho, beta];

% --- 2. Benettin's Algorithm Parameters ---
% Initial conditions for the fiducial trajectory
X_0 = [-4.7795941536989090537     -8.8325236589201416848      10.425932436960884786]'; % X_0

% Initial small perturbation
d0 = 1e-4; % Initial separation magnitude


% --- 3. Time and sampling parameters
% Time parameters for saving and plotting the state
fs = 100; % plotting sample frequency
dt = 1/fs; % plotting sample time
T = [-100, 200]; % plotting time interval, from T(1) to 0, is warmup and settle down which isn't typically plotted.
if not( T(1)<=0 & 0<T(2) )
    error('T(1) must be 0 or negative, and T(2) must be positive')
end
nt = (T(2)-T(1))*fs+1; % number of plotting samples
t = linspace(T(1), T(2), nt)'; % plotting time vector

lya_dt = 0.1; % seconds

%% --- 4. ODE solver options and system integration
ode_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13, 'MaxStep', 0.1*dt, 'InitialStep', 0.05*dt);   % never step more than 10 % of dt

% Integrate the system along t with X_0 ICs
[t_ode, X] = ode45(@(tt,XX) lorenz_eqs(tt,XX,params), t, X_0, ode_options);
assert(all(abs(t_ode - t) < 1e-14), 'ODE solver did not return results exactly at the requested times.');
clear t_ode

%% --- 5. Benettin's algorithm call -----------------
[LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t, dt, fs, d0, T, lya_dt, params, ode_options, @lorenz_eqs);

fprintf('----------------------------------------------------\n');
fprintf('Estimated Largest Lyapunov Exponent (LLE): %f\n', LLE);
fprintf('----------------------------------------------------\n');

% --- Plot State Variables and Lyapunov Exponents ---
figure;
ax = gobjects(4,1); % Preallocate array for axes handles

% Subplot 1: Time vs. x-state
ax(1) = subplot(4,1,1);
plot(t, X(:,1));
xlabel('Time (t)');
ylabel('x(t)');
title('State Variable x');
grid on;

% Subplot 2: Time vs. y-state
ax(2) = subplot(4,1,2);
plot(t, X(:,2));
xlabel('Time (t)');
ylabel('y(t)');
title('State Variable y');
grid on;

% Subplot 3: Time vs. z-state
ax(3) = subplot(4,1,3);
plot(t, X(:,3));
xlabel('Time (t)');
ylabel('z(t)');
title('State Variable z');
grid on;

% Subplot 4: Time vs. Lyapunov Exponents
ax(4) = subplot(4,1,4);
plot(t_lya, local_lya, 'DisplayName', 'Local LLE');
hold on;
plot(t_lya, finite_lya, 'DisplayName', 'Finite-time LLE');
plot([t_lya(1) t_lya(end)], [LLE LLE], 'r--', 'DisplayName', sprintf('Final LLE: %.4f', LLE));
hold off;
xlabel('Time (t)');
ylabel('Lyapunov Exponent');
title('Lyapunov Exponent Estimates');
legend('show', 'Location', 'best');
grid on;

linkaxes(ax, 'x'); % Link x-axes of all subplots


% --- 6. Plot convergence of LLE (optional) ---
figure;
plot(t_lya, finite_lya);
xlabel('Total Time (t)');
ylabel('Largest Lyapunov Exponent Estimate (\lambda_1)');
title('Convergence of the Largest Lyapunov Exponent');
grid on;
hold on;
plot([0, T(2)], [LLE, LLE], 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Final LLE: %.4f', LLE));
legend('show', 'Location', 'best');
hold off;

% --- Lorenz System Equations Function ---
function dXdt = lorenz_eqs(~, X, params)
    % Lorenz system differential equations
    % X = [x; y; z] state vector
    % params = [sigma; rho; beta] parameters

    sigma = params(1);
    rho   = params(2);
    beta  = params(3);

    x = X(1);
    y = X(2);
    z = X(3);

    dXdt = zeros(3,1);
    dXdt(1) = sigma * (y - x);
    dXdt(2) = x * (rho - z) - y;
    dXdt(3) = x * y - beta * z;
end

% --- Benettin's Algorithm ---
function [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t, dt, fs, d0, T, lya_dt, params, ode_options, dynamics_func)
    % reshoots small segments to compute the divergence rate along the system trajectory in X
    deci_lya = round(lya_dt*fs);     % samples per Lyapunov interval
    tau_lya = dt*deci_lya;    % Lyapunov rescaling time interval (integration time between rescalings)
    t_lya    = t(1:deci_lya:end);     % direct decimation of `t`

    if t_lya(end) + tau_lya > T(2)    % keep segments fully inside [T(1),T(2)]
        t_lya(end) = [];
    end
    nt_lya  = numel(t_lya);           % number of Lyapunov intervals

    local_lya  = zeros(nt_lya,1);
    finite_lya = zeros(nt_lya,1);
    sum_log_stretching_factors = 0;

    % Initial perturbation
    rnd_IC = randn(3,1);
    pert = (rnd_IC./norm(rnd_IC)).*d0;

    for k = 1:nt_lya
        idx_start = (k-1)*deci_lya + 1;      % index of t_lya(k) in `t`
        idx_end   = idx_start + deci_lya;    % => t_lya(k) + tau_lya

        X_start = X(idx_start,:).';      % fiducial state at t_lya(k)

        X_k_pert = X_start + pert; % Rescale perturbation from previous normalized delta

        % Integrate ONLY the perturbed trajectory over [t_lya(k), t_lya(k)+tau_lya]
        t_seg = t_lya(k) + [0, tau_lya];
        [~, X_pert_seg] = ode45(@(tt,XX) dynamics_func(tt,XX,params), t_seg, X_k_pert, ode_options);
        
        X_pert_end = X_pert_seg(end,:).';

        X_end   = X(idx_end,:).';        % fiducial state tau_lya later

        % Local exponent
        delta   = X_pert_end - X_end;
        d_k     = norm(delta);
        local_lya(k) = log(d_k/d0)/tau_lya;

        pert = (delta./d_k).*d0; % rescaled perturbation for next step

        % Accumulate stretching only from t >= 0
        if t_lya(k) >= 0
            sum_log_stretching_factors = sum_log_stretching_factors + log(d_k/d0);
            finite_lya(k,1) = sum_log_stretching_factors / max(t_lya(k)+tau_lya, eps);
        end
    end

    LLE = sum_log_stretching_factors / T(2);  % finite-time estimate from t = 0 to T(2)
end