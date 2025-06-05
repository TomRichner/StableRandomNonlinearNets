function SRNN_tseries_plot(t, u_ex, r, a, b, u_d, params, EI_vec, T, n, M, Lya_method, varargin)
% SRNN_tseries_plot - Plot time series results from SRNN simulation
%
% Inputs:
%   t        - Time vector
%   u_ex     - External input (n x nt)
%   r        - Firing rates (n x nt)
%   a        - SFA variables ((n*n_a) x nt)
%   b        - STD variables ((n*n_b) x nt) 
%   u_d      - Dendritic potential (n x nt)
%   params   - Parameter structure containing n_a, n_b, c_SFA, F_STD
%   EI_vec   - Excitatory/inhibitory classification vector (n x 1)
%   T        - Time interval [T_start, T_end]
%   n        - Number of neurons
%   M        - Connectivity matrix (n x n)
%   Lya_method - Lyapunov calculation method ('benettin', 'qr', or 'none')
%   varargin - Optional Lyapunov results: LLE, t_lya, local_lya, finite_lya, 
%              LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, N_sys_eqs

% Parse optional Lyapunov inputs
LLE = [];
t_lya = [];
local_lya = [];
finite_lya = [];
LE_spectrum = [];
local_LE_spectrum_t = [];
finite_LE_spectrum_t = [];
N_sys_eqs = [];

if ~isempty(varargin)
    lya_results = varargin{1};
    if isfield(lya_results, 'LLE'), LLE = lya_results.LLE; end
    if isfield(lya_results, 't_lya'), t_lya = lya_results.t_lya; end
    if isfield(lya_results, 'local_lya'), local_lya = lya_results.local_lya; end
    if isfield(lya_results, 'finite_lya'), finite_lya = lya_results.finite_lya; end
    if isfield(lya_results, 'LE_spectrum'), LE_spectrum = lya_results.LE_spectrum; end
    if isfield(lya_results, 'local_LE_spectrum_t'), local_LE_spectrum_t = lya_results.local_LE_spectrum_t; end
    if isfield(lya_results, 'finite_LE_spectrum_t'), finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t; end
    if isfield(lya_results, 'N_sys_eqs'), N_sys_eqs = lya_results.N_sys_eqs; end
end

nt = length(t);

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
ylabel('Input, u_{ex}');
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
ylabel('Spike rate, r');
box off;
set(gca, 'XTickLabel', []);

% Subplot 3: SFA sum (from variable 'a')
ax_handles(3) = subplot(num_subplots, 1, 3);
if params.n_a > 0
    % a is n x n_a x nt from unpack_SRNN_state when called with full X
    % sum(a, 2) results in n x 1 x nt. We want n x nt.
    a_sum_plot = squeeze_dim(sum(a, 2), 2); % Use squeeze_dim to remove singleton 2nd dim
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
ylabel({'Spike Freq. Adapt.', '$\sum\limits_k a_k$'}, 'Interpreter', 'latex')
box off;
set(gca, 'XTickLabel', []);

% Subplot 4: STD product (from variable 'b')
ax_handles(4) = subplot(num_subplots, 1, 4);
if params.n_b > 0
    % b is n x n_b x nt from unpack_SRNN_state when called with full X
    % prod(b, 2) results in n x 1 x nt. We want n x nt.
    b_prod_plot = squeeze_dim(prod(b, 2), 2); % Use squeeze_dim to remove singleton 2nd dim
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
ylabel({'Syn. Dep.','$\prod\limits_m b_m$'}, 'Interpreter', 'latex')
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
ylabel('Dentrite, u_d');
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
            plot([t_lya(1) t_lya(end)], [LLE LLE], 'Color',[0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', sprintf('Global LLE: %.4f', LLE));
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

% Network graph visualization (if network is small enough)
if n <= 50
    figure(2)
    clf
    set(gcf,'Position',[100   1009 500 500])
    [h_digraph, dgA] = plot_network_graph_widthRange_color_R(M,1,EI_vec);
    box off
    axis equal
    axis off
end

end 