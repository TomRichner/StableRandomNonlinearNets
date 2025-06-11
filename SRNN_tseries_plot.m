function SRNN_tseries_plot(t, u_ex, r_ts, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T_plot_limits, Lya_method, varargin)
% SRNN_tseries_plot - Plot time series results from SRNN simulation
%
% Inputs:
%   t        - Time vector
%   u_ex     - External input (n x nt)
%   r_ts     - Firing rates (n x nt)
%   a_E_ts   - SFA variables for E neurons (n_E x n_a_E x nt or empty)
%   a_I_ts   - SFA variables for I neurons (n_I x n_a_I x nt or empty)
%   b_E_ts   - STD variables for E neurons (n_E x n_b_E x nt or empty)
%   b_I_ts   - STD variables for I neurons (n_I x n_b_I x nt or empty)
%   u_d_ts   - Dendritic potential (n x nt)
%   params   - Parameter structure
%   T_plot_limits - Time interval [T_start, T_end] for x-axis limits
%   Lya_method - Lyapunov calculation method ('benettin', 'qr', or 'none')
%   varargin - Optional Lyapunov results struct

% Parse optional Lyapunov inputs
LLE = [];
t_lya = [];
local_lya = [];
finite_lya = [];
LE_spectrum = [];
local_LE_spectrum_t = [];
finite_LE_spectrum_t = [];
N_sys_eqs_lya = []; % Renamed to avoid conflict with params.N_sys_eqs if it exists

if ~isempty(varargin)
    lya_results = varargin{1};
    if isfield(lya_results, 'LLE'), LLE = lya_results.LLE; end
    if isfield(lya_results, 't_lya'), t_lya = lya_results.t_lya; end
    if isfield(lya_results, 'local_lya'), local_lya = lya_results.local_lya; end
    if isfield(lya_results, 'finite_lya'), finite_lya = lya_results.finite_lya; end
    if isfield(lya_results, 'LE_spectrum'), LE_spectrum = lya_results.LE_spectrum; end
    if isfield(lya_results, 'local_LE_spectrum_t'), local_LE_spectrum_t = lya_results.local_LE_spectrum_t; end
    if isfield(lya_results, 'finite_LE_spectrum_t'), finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t; end
    if isfield(lya_results, 'N_sys_eqs'), N_sys_eqs_lya = lya_results.N_sys_eqs; end % Use from lya_results
end

nt = length(t);
n = params.n;
E_indices = params.E_indices; % from params struct
I_indices = params.I_indices; % from params struct
M = params.M; % from params struct


%% Make plots

figure('Position', [100, 100, 900, 1200]);

num_subplots = 5;
if ~strcmpi(Lya_method, 'none')
    num_subplots = 6;
end
ax_handles = gobjects(num_subplots, 1);

plot_indices = round(linspace(1, nt, min(nt, 2000)));
t_display = t(plot_indices);

% Subplot 1: External input u_ex
ax_handles(1) = subplot(num_subplots, 1, 1);
has_stim = any(abs(u_ex(:, plot_indices)) > 1e-6, 2); 
if any(has_stim)
    plot(t_display, u_ex(has_stim, plot_indices)');
else
    plot(t_display, zeros(length(t_display),1)); 
end
ylabel('Input, u_{ex}');
box off;
set(gca, 'XTickLabel', []);

% Subplot 2: Firing rates r_ts
ax_handles(2) = subplot(num_subplots, 1, 2);
hold on;
if ~isempty(I_indices)
    plot(t_display, r_ts(I_indices, plot_indices)', 'r');
end
if ~isempty(E_indices)
    set(gca,'ColorOrderIndex',1); 
    plot(t_display, r_ts(E_indices, plot_indices)');
end
hold off;
ylabel('Spike rate, r');
box off;
set(gca, 'XTickLabel', []);

% Subplot 3: SFA sum
ax_handles(3) = subplot(num_subplots, 1, 3);
a_sum_plot = zeros(n, nt); % Initialize as n x nt
if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
    % a_E_ts is n_E x n_a_E x nt. sum across 2nd dim -> n_E x 1 x nt. Squeeze -> n_E x nt
    a_sum_plot(E_indices, :) = squeeze(sum(a_E_ts, 2));
end
if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
    a_sum_plot(I_indices, :) = squeeze(sum(a_I_ts, 2));
end

hold on;
% Plot I neurons SFA sum if they have SFA states and c_SFA is non-zero for them
if ~isempty(I_indices) && params.n_a_I > 0
    active_I_sfa = params.c_SFA(I_indices) ~= 0;
    if any(active_I_sfa)
         plot(t_display, a_sum_plot(I_indices(active_I_sfa), plot_indices)', 'r');
    end
end
% Plot E neurons SFA sum
if ~isempty(E_indices) && params.n_a_E > 0
    active_E_sfa = params.c_SFA(E_indices) ~= 0;
    if any(active_E_sfa)
        set(gca,'ColorOrderIndex',1);
        plot(t_display, a_sum_plot(E_indices(active_E_sfa), plot_indices)');
    end
end
hold off;
ylabel({'Spike Freq. Adapt.', '$\sum\limits_k a_k$'}, 'Interpreter', 'latex')
box off;
set(gca, 'XTickLabel', []);

% Subplot 4: STD product
ax_handles(4) = subplot(num_subplots, 1, 4);
b_prod_plot = ones(n, nt); % Initialize as n x nt, default product is 1
if params.n_E > 0 && params.n_b_E > 0 && ~isempty(b_E_ts)
    % b_E_ts is n_E x n_b_E x nt. prod across 2nd dim -> n_E x 1 x nt. Squeeze -> n_E x nt
    b_prod_plot(E_indices, :) = squeeze(prod(b_E_ts, 2));
end
if params.n_I > 0 && params.n_b_I > 0 && ~isempty(b_I_ts)
    b_prod_plot(I_indices, :) = squeeze(prod(b_I_ts, 2));
end

hold on;
% Plot I neurons STD product if they have STD states and F_STD is non-zero
if ~isempty(I_indices) && params.n_b_I > 0
    active_I_std = params.F_STD(I_indices) ~= 0;
    if any(active_I_std)
        plot(t_display, b_prod_plot(I_indices(active_I_std), plot_indices)', 'r');
    end
end
% Plot E neurons STD product
if ~isempty(E_indices) && params.n_b_E > 0
    active_E_std = params.F_STD(E_indices) ~= 0;
    if any(active_E_std)
        set(gca,'ColorOrderIndex',1);
        plot(t_display, b_prod_plot(E_indices(active_E_std), plot_indices)');
    end
end
hold off;
ylabel({'Syn. Dep.','$\prod\limits_m b_m$'}, 'Interpreter', 'latex')
box off;
ylim([0 1.1]); 
set(gca, 'XTickLabel', []);

% Subplot 5: Dendritic potential u_d_ts
ax_handles(5) = subplot(num_subplots, 1, 5);
hold on;
if ~isempty(I_indices)
    plot(t_display, u_d_ts(I_indices, plot_indices)', 'r');
end
if ~isempty(E_indices)
    set(gca,'ColorOrderIndex',1);
    plot(t_display, u_d_ts(E_indices, plot_indices)');
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
        % Check if there is any Lyapunov data to plot.
        % t_lya would be empty if divergence occurred before the first interval.
        if ~isempty(t_lya)
            plot([t_lya(1) t_lya(end)], [LLE LLE], 'Color',[0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', sprintf('Global LLE: %.4f', LLE));
            plot(t_lya, finite_lya, 'r', 'LineWidth', 2, 'DisplayName', 'Finite-time LLE');
            plot(t_lya, local_lya, 'b', 'LineWidth', 1, 'DisplayName', 'Local LLE');
            
            % Adjust y-limits to ensure all data is visible
            ylim_current = get(gca, 'YLim');
            ylim_new = [min(ylim_current(1), LLE - 0.05*abs(LLE)), max(ylim_current(2), LLE + 0.05*abs(LLE))];
            if all(isfinite(ylim_new))
                ylim(ylim_new);
            end
        else
            text(0.5, 0.5, 'No Lyapunov data (diverged early)', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end

    elseif strcmpi(Lya_method, 'qr')
        if ~isempty(LE_spectrum) && ~isempty(t_lya) && ~isempty(N_sys_eqs_lya)
            colors = lines(N_sys_eqs_lya);
            if ~isempty(local_LE_spectrum_t)
                for i = 1:N_sys_eqs_lya
                    plot(t_lya, local_LE_spectrum_t(:,i), '--', 'Color', colors(i,:), 'DisplayName', sprintf('Local LE(%d)', i));
                end
            end
            if ~isempty(finite_LE_spectrum_t)
                for i = 1:N_sys_eqs_lya
                    plot(t_lya, finite_LE_spectrum_t(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Finite LE(%d)', i));
                end
            end
            for i = 1:N_sys_eqs_lya
                plot([t_lya(1) t_lya(end)], [LE_spectrum(i) LE_spectrum(i)], ':', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('Global LE(%d): %.4f', i, LE_spectrum(i)));
            end
            % For QR method, fix y-limits for better comparability of spectra
            ylim([-0.5 0.5]);
        else
             text(0.5, 0.5, 'No Lyapunov data to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
    end
    hold off;
    ylabel('Lyapunov Exp.');
    xlabel('Time (s)');
    legend('show', 'Location', 'best');
    grid on;
    box off;
    ylim([-0.5 0.5])
end

% Link all subplots' x-axes and set limits
if numel(ax_handles) > 0 && all(isgraphics(ax_handles))
    linkaxes(ax_handles, 'x');
    % Set xlim on the first axes, and linkaxes will propagate it to the others.
    axes(ax_handles(1)); 
    xlim(T_plot_limits);
end

if n <= 50
    figure(2)
    clf
    set(gcf,'Position',[100   1009 500 500])
    % EI_vec is already in params, so params.EI_vec
    [h_digraph, dgA] = plot_network_graph_widthRange_color_R(M,1,params.EI_vec);
    box off
    axis equal
    axis off
end

end 