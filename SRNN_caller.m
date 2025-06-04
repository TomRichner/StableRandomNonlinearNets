% example call of SRRN()

close all
clear
clc

tic

%% 
seed = 7;
rng(seed,'twister');

%% Network
n = 10; % number of neurons

Lya_method = 'benettin'; % 'benettin', 'qr', or 'none

mean_in_out_degree = 4; % desired mean number of connections in and out
density = mean_in_out_degree/(n-1); % each neuron can make up to n-1 connections with other neurons
sparsity = 1-density;

EI = 0.7;
scale = 0.5/0.79782; % overall scaling factor of weights
w.EE = scale*2; % E to E. Change to scale*2 for bursting
w.EI = scale*1; % E to I connections
w.IE = scale*1; % I to E
w.II = scale*.5; % I to I
w.selfE = 0;    % self connections of E neurons
w.selfI = 0;    % self connections of I neurons

[M, EI_vec] = generate_M(n,w,sparsity, EI);
EI_vec = EI_vec(:); % make it a column

%% Time
fs = 250; %Plotting sample frequency
dt = 1/fs;
T = [-10 100];

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
u_ex(1,-t(1)*fs+fs*1+(1:fix(fs*dur))) = stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));
u_ex(1,-t(1)*fs+fix(fs*5)+(1:fix(fs*dur))) = stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))');
u_ex = u_ex*1;
u_ex = u_ex(:,1:nt);
DC = 0.001;
u_ex = u_ex+DC;
% u_ex(:,1:fs) = u_ex(:,1:fs)+10./fs.*randn(n,fs); % noise in the first second to help the network get off the trivial saddle node from ICs
u_ex = u_ex+0.001./fs.*randn(n,nt); % a tiny bit of noise to help the network get off the trivial saddle node from ICs


%% parameters

n_a = 3; % number of SFA timescales per neuron
n_b = 1; % number of STD timescales per neuron

tau_a = logspace(log10(0.3), log10(6), n_a); % s, 1 x n_a, time constants of SFA
tau_b = logspace(log10(0.6), log10(9), n_b);  % s, 1 x n_b, time constants of STD, n_b == 1, then it takes the last value log10(9)

tau_d = 0.01; % s, scalar

c_SFA = 1 * double(EI_vec == 1); % n x 1, 0 for I neurons bc no SFA
F_STD = 1 * double(EI_vec == 1); % n x 1, 0 for I neurons bc no STD
tau_STD = 1; % scalar, time constant of synaptic depression

params = package_params(n_a, n_b, tau_a, tau_b, tau_d, n, M, c_SFA, F_STD, tau_STD);

%% Initial Conditions
if n_a > 0
    a0 = zeros(n * n_a, 1);
else
    a0 = [];
end

if n_b > 0
    b0 = ones(n * n_b, 1);
else
    b0 = [];
end

u_d0 = zeros(n, 1);

X_0 = [a0; b0; u_d0];

N_sys_eqs = size(X_0,1); % Number of system equations / states

%% Integrate with ODE solver

ode_options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt, 'InitialStep', min(0.001, 0.2*dt));

SRNN_wrapper = @(tt,XX) SRNN(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params

[t_ode, X] = ode15s(SRNN_wrapper, t, X_0, ode_options);

assert(all(abs(t_ode - t) < 1e-12), 'ODE solver did not return results exactly at the requested times for fiducial trajectory.');
clear t_ode % t_ode is same as t

%% comput LLE or Lyapunov spectrum

lya_dt = 0.1; % Rescaling time interval for Lyapunov calculation (tau_lya) (s)

switch lower(Lya_method)
    case 'qr'
        fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
        
        % Ensure SRNN_jacobian_eqs is defined elsewhere or this will error.
        % Using N_sys_eqs for the number of states.
        [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = ...
            lyapunov_spectrum_qr(X, t, lya_dt, params, ode_options, @SRNN_jacobian_eqs, T, N_sys_eqs, fs);

        % Display the estimated Lyapunov Spectrum
        fprintf('----------------------------------------------------\n');
        fprintf('Estimated Lyapunov Spectrum (Global):\n');
        for i = 1:N_sys_eqs
            fprintf('  LE(%d): %f\n', i, LE_spectrum(i));
        end
        fprintf('Sum of exponents: %f (should be < 0 for dissipative systems like Lorenz)\n', sum(LE_spectrum));
        fprintf('Kaplan-Yorke Dimension: %f\n', calculate_kaplan_yorke_dimension(LE_spectrum));
        fprintf('----------------------------------------------------\n');
        
    case 'benettin'
        fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
        
        d0 = 1e-4; % Initial separation magnitude for Benettin's algorithm
        [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex);

        fprintf('----------------------------------------------------\n');
        fprintf('Estimated Largest Lyapunov Exponent (LLE): %f\n', LLE);
        fprintf('----------------------------------------------------\n');
        
    case 'none'
        fprintf('Skipping Lyapunov calculation - trajectory only.\n');
        
    otherwise
        error('Unknown method: %s. Choose ''qr'', ''benettin'', or ''none''.', method);
end

%% Convert X to named variables
[a, b, u_d] = unpack_SRNN_state(X, n, n_a, n_b);

% compute dependent variables r and p
%% Compute dependent variables r and p using a subfunction
[r, p] = compute_dependent_variables(a, b, u_d, n_a, n_b, c_SFA);

%% Make plots

figure('Position', [100, 100, 900, 1200]); % Adjusted figure size

num_subplots = 5;
if ~strcmpi(Lya_method, 'none')
    num_subplots = 6;
end
ax_handles = gobjects(num_subplots, 1);

% Define plot indices for manageable number of points
plot_indices = round(linspace(1, nt, min(nt, 2000)));
t_display = t(plot_indices);

E_neurons_idx = find(EI_vec == 1);
I_neurons_idx = find(EI_vec == -1);

% Subplot 1: External input u_ex
ax_handles(1) = subplot(num_subplots, 1, 1);
has_stim = any(abs(u_ex) > 1e-6, 2); % Find neurons that actually receive stimulus
if any(has_stim)
    plot(t_display, u_ex(has_stim, plot_indices)');
else
    plot(t_display, zeros(length(t_display),1)); % Plot zeros if no stim
end
ylabel('u_{ex}');
box off;
set(gca, 'XTickLabel', []);

% Subplot 2: Firing rates r
ax_handles(2) = subplot(num_subplots, 1, 2);
hold on;
if ~isempty(I_neurons_idx)
    plot(t_display, r(I_neurons_idx, plot_indices)', 'r');
end
if ~isempty(E_neurons_idx)
    set(gca,'ColorOrderIndex',1); % Reset color index if I neurons were plotted
    plot(t_display, r(E_neurons_idx, plot_indices)');
end
hold off;
ylabel('r');
box off;
set(gca, 'XTickLabel', []);

% Subplot 3: SFA sum (from variable 'a')
ax_handles(3) = subplot(num_subplots, 1, 3);
if params.n_a > 0
    a_reshaped_plotting = reshape(a, n, params.n_a, nt); % a is (n*n_a) x nt
    a_sum_plot = squeeze(sum(a_reshaped_plotting, 2)); % n x nt
else
    a_sum_plot = zeros(n, nt);
end
hold on;
if ~isempty(I_neurons_idx) && params.n_a > 0 % Assuming SFA can apply to I neurons based on c_SFA
    if any(params.c_SFA(I_neurons_idx) ~= 0) % Only plot if SFA is active for I neurons
         plot(t_display, a_sum_plot(I_neurons_idx, plot_indices)', 'r');
    end
end
if ~isempty(E_neurons_idx) && params.n_a > 0
    set(gca,'ColorOrderIndex',1);
    plot(t_display, a_sum_plot(E_neurons_idx, plot_indices)');
end
hold off;
ylabel('$\\sum_k a_k$ (SFA)', 'Interpreter', 'latex');
box off;
set(gca, 'XTickLabel', []);

% Subplot 4: STD product (from variable 'b')
ax_handles(4) = subplot(num_subplots, 1, 4);
if params.n_b > 0
    b_reshaped_plotting = reshape(b, n, params.n_b, nt); % b is (n*n_b) x nt
    b_prod_plot = squeeze(prod(b_reshaped_plotting, 2)); % n x nt
else
    b_prod_plot = ones(n, nt);
end
hold on;
if ~isempty(I_neurons_idx) && params.n_b > 0 % Assuming STD can apply to I neurons based on F_STD
    if any(params.F_STD(I_neurons_idx) ~= 0) % Only plot if STD is active for I neurons
        plot(t_display, b_prod_plot(I_neurons_idx, plot_indices)', 'r');
    end
end
if ~isempty(E_neurons_idx) && params.n_b > 0
    set(gca,'ColorOrderIndex',1);
    plot(t_display, b_prod_plot(E_neurons_idx, plot_indices)');
end
hold off;
ylabel('$\\prod_k b_k$ (STD)', 'Interpreter', 'latex');
box off;
ylim([0 1.1]); % STD factors are typically between 0 and 1
set(gca, 'XTickLabel', []);

% Subplot 5: Dendritic potential u_d
ax_handles(5) = subplot(num_subplots, 1, 5);
hold on;
if ~isempty(I_neurons_idx)
    plot(t_display, u_d(I_neurons_idx, plot_indices)', 'r');
end
if ~isempty(E_neurons_idx)
    set(gca,'ColorOrderIndex',1);
    plot(t_display, u_d(E_neurons_idx, plot_indices)');
end
hold off;
ylabel('u_d');
box off;
if num_subplots == 5
    xlabel('Time (s)');
else
    set(gca, 'XTickLabel', []);
end

% Subplot 6: Lyapunov Exponents (if calculated)
if ~strcmpi(Lya_method, 'none')
    ax_handles(6) = subplot(num_subplots, 1, 6);
    hold on;
    if strcmpi(Lya_method, 'benettin')
        if exist('LLE', 'var') && exist('t_lya', 'var') && ~isempty(t_lya)
            plot([t_lya(1) t_lya(end)], [LLE LLE], 'k', 'LineWidth', 3, 'DisplayName', sprintf('Global LLE: %.4f', LLE));
        end
            if exist('t_lya', 'var') && exist('finite_lya', 'var') && ~isempty(t_lya)
            plot(t_lya, finite_lya, 'r', 'LineWidth', 2, 'DisplayName', 'Finite-time LLE');
        end
        if exist('t_lya', 'var') && exist('local_lya', 'var') && ~isempty(t_lya)
            plot(t_lya, local_lya, 'b', 'LineWidth', 1, 'DisplayName', 'Local LLE');
        end
        ylim([-.5 .5])
    elseif strcmpi(Lya_method, 'qr')
        if exist('LE_spectrum', 'var') && exist('t_lya', 'var') && ~isempty(t_lya)
            colors = lines(N_sys_eqs);
            % Plot local Lyapunov exponents
            if exist('local_LE_spectrum_t', 'var') && ~isempty(local_LE_spectrum_t)
                for i = 1:N_sys_eqs
                    plot(t_lya, local_LE_spectrum_t(:,i), '--', 'Color', colors(i,:), 'DisplayName', sprintf('Local LE(%d)', i));
                end
            end
            % Plot finite-time Lyapunov exponents
            if exist('finite_LE_spectrum_t', 'var') && ~isempty(finite_LE_spectrum_t)
                for i = 1:N_sys_eqs
                    plot(t_lya, finite_LE_spectrum_t(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Finite LE(%d)', i));
                end
            end
            % Plot global Lyapunov exponents as horizontal lines
            for i = 1:N_sys_eqs
                plot([t_lya(1) t_lya(end)], [LE_spectrum(i) LE_spectrum(i)], ':', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('Global LE(%d): %.4f', i, LE_spectrum(i)));
            end
        else
             plot(0,0); % Placeholder if Lya vars are missing
        end
    end
    hold off;
    ylabel('Lyapunov Exp.');
    xlabel('Time (s)');
    legend('show');
    grid on;
    box off;
end

% Link axes and set limits
linkaxes(ax_handles, 'x');
xlim(T);

if n <= 50
    figure(2)
    clf
    set(gcf,'Position',[100   1009 500 500])
    [h_digraph, dgA] = plot_network_graph_widthRange_color_R(M,1,EI_vec);
    box off
    axis equal
    axis off
end

sim_dur = toc

sim_t_dived_by_rt = sim_dur./(T(2)-T(1))