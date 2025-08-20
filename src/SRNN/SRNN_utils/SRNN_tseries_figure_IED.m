function SRNN_tseries_figure_IED(t, u_ex, r_ts, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T_plot_limits, Lya_method, sr_or_poisson, include_sum_E_I_SR, varargin)
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
%   sr_or_poisson - 'sr' for spike rate plot, 'poisson' for raster plot
%   include_sum_E_I_SR - boolean to control plotting of sum E/I spike rates
%   varargin - Optional Lyapunov results struct

% --- Plotting options ---

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

% Set default for sr_or_poisson if not provided
if nargin < 12 || isempty(sr_or_poisson)
    sr_or_poisson = 'sr';
end

scale_bar_info = []; % Store scale bar parameters here

%% Make plots

figure('Position', [100, 100, 870, 1130]);

% Calculate total number of subplots
num_subplots_base = 4; % Input, Rates, SFA, STD
if include_sum_E_I_SR
    num_subplots_base = num_subplots_base + 1;
end
if ~strcmpi(Lya_method, 'none')
    num_subplots_base = num_subplots_base + 1;
end

plot_offset = 0;
if strcmpi(sr_or_poisson, 'poisson') || strcmpi(sr_or_poisson, 'sr_stacked')
    plot_offset = 1;
end
num_subplots = num_subplots_base + plot_offset;

ax_handles = gobjects(0);
sp_idx = 1;

plot_indices = round(linspace(1, nt, min(nt, 2000)));
t_display = t(plot_indices);

% Subplot 1: External input u_ex
ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
sp_idx = sp_idx + 1;
has_stim = any(abs(u_ex(:, plot_indices)) > 1e-6, 2); 
if any(has_stim)
    hold on;
    
    inh_color_rgb = [1 0 0]; % Red for inhibitory

    % Plot inhibitory neurons with stimulus
    if ~isempty(I_indices)
        for k = 1:numel(I_indices)
            i = I_indices(k);
            if has_stim(i)
                plot(t_display, u_ex(i, plot_indices)', 'Color', inh_color_rgb);
            end
        end
    end

    % Plot excitatory neurons with stimulus
    if ~isempty(E_indices)
        cmap_exc = lines(numel(E_indices));
        for k = numel(E_indices):-1:1
            i = E_indices(k);
            if has_stim(i)
                plot_color = cmap_exc(k, :);
                plot(t_display, u_ex(i, plot_indices)', 'Color', plot_color);
            end
        end
    end

    hold off;
else
    plot(t_display, zeros(length(t_display),1)); 
end
ylabel('Input');
box off;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

% Subplot 2: Sum E and I spike rates
if include_sum_E_I_SR
    ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
    sp_idx = sp_idx + 1;
    hold on;
    if ~isempty(E_indices)
        % sum_E_rate = sum(r_ts(E_indices, plot_indices), 1);
        % plot(t_display, sum_E_rate, 'b', 'LineWidth', 2, 'DisplayName', 'excitatory neurons, sum rate');
        mean_E_rate = mean(r_ts(E_indices, plot_indices), 1);
        plot(t_display, mean_E_rate, 'b', 'LineWidth', 2, 'DisplayName', 'excitatory neurons');
    end
    if ~isempty(I_indices)
        % sum_I_rate = sum(r_ts(I_indices, plot_indices), 1);
        % plot(t_display, sum_I_rate, 'r', 'LineWidth', 2, 'DisplayName', 'inhibitory neurons');
        mean_I_rate = mean(r_ts(I_indices, plot_indices), 1);
        plot(t_display, mean_I_rate, 'r', 'LineWidth', 2, 'DisplayName', 'inhibitory neurons');
    end
    hold off;
    ylabel('Mean Firing Rate');
    box off;
    set(gca, 'XTickLabel', []);
    set(gca, 'XTick', []);
    leg = legend;
    leg_pos = leg.Position;
    % leg.Position = [leg_pos(1) - 0.4, leg_pos(2) + 0.031, leg_pos(3), leg_pos(4)]; % move this legend up and left a bit
    leg.Position = [leg_pos(1) - 0.0, leg_pos(2) + 0.031, leg_pos(3), leg_pos(4)]; % move this legend up and left a bit
    leg.Box = 'off';
end

% Subplot 3: Firing rates r_ts or Poisson raster plot
if strcmpi(sr_or_poisson, 'poisson') || strcmpi(sr_or_poisson, 'sr_stacked')
    ax_handles(end+1) = subplot(num_subplots, 1, sp_idx:(sp_idx + plot_offset));
    sp_idx = sp_idx + 1 + plot_offset;
else
    ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
    sp_idx = sp_idx + 1;
end

if strcmpi(sr_or_poisson, 'poisson')
    % Poisson raster plot implementasbtion
    step_factor = 2;  % sample every dt for Poisson draw
    dt = t(2) - t(1);  % time step from input
    dt_poisson = step_factor * dt;
    
    % Down-sample the spike-rate matrix
    idx_sample = 1:step_factor:nt;
    t_sample = t(idx_sample);
    r_sample = r_ts(:, idx_sample);
    
    % Poisson draw: number of spikes per bin ~ Poisson(lambda)
    lambda = r_sample * dt_poisson;  % expected spikes per bin
    spike_counts = poissrnd(lambda);
    
    % Find all (neuron, time) bins with ≥1 spike
    [spk_neuron, spk_time_idx] = find(spike_counts > 0);
    spk_times = t_sample(spk_time_idx);
    counts_in_bin = spike_counts(sub2ind(size(spike_counts), spk_neuron, spk_time_idx));
    
    % Dither extra spikes in time instead of making taller ticks
    % Create arrays to hold all individual spike times and neuron IDs
    all_spike_times = [];
    all_spike_neurons = [];
    
    dither_amount = dt_poisson * 10; % Dither within ±10x of the bin width
    
    for i = 1:length(spk_times)
        n_spikes = counts_in_bin(i);
        base_time = spk_times(i);
        neuron_id = spk_neuron(i);
        
        if n_spikes == 1
            % Single spike - no dithering needed
            all_spike_times = [all_spike_times; base_time];
            all_spike_neurons = [all_spike_neurons; neuron_id];
        else
            % Multiple spikes - dither them around the base time
            dithered_times = base_time + dither_amount * randn(n_spikes, 1);
            all_spike_times = [all_spike_times; dithered_times];
            all_spike_neurons = [all_spike_neurons; repmat(neuron_id, n_spikes, 1)];
        end
    end
    
    % Plot parameters for vertical ticks
    tick_base_height_val = 0.75;  % Base height for a single spike
    tick_line_width = 2;
    inh_color_rgb = [1 0 0];  % Red for inhibitory
    
    % Get colormap for excitatory neurons
    if ~isempty(E_indices)
        cmap_exc = lines(numel(E_indices));
    else
        cmap_exc = [];
    end
    
    hold on;
    
    % Plot inhibitory neurons first (red)
    if ~isempty(I_indices)
        for k_loop_inh = 1:numel(I_indices)
            n_neuron = I_indices(k_loop_inh);
            spike_indices_for_neuron = find(all_spike_neurons == n_neuron);
            
            for i_spike = 1:length(spike_indices_for_neuron)
                spike_idx = spike_indices_for_neuron(i_spike);
                x_time = all_spike_times(spike_idx);
                y_center = n_neuron;
                
                % All spikes have the same height now
                current_tick_h = tick_base_height_val;
                y_coords = [y_center - current_tick_h/2, y_center + current_tick_h/2];
                
                plot([x_time x_time], y_coords, ...
                     'Color', inh_color_rgb, ...
                     'LineWidth', tick_line_width);
            end
        end
    end
    
    % Plot excitatory neurons with lines colormap
    if ~isempty(E_indices)
        for k_loop_exc = 1:numel(E_indices)
            n_neuron = E_indices(k_loop_exc);
            spike_indices_for_neuron = find(all_spike_neurons == n_neuron);
            neuron_color_rgb = cmap_exc(k_loop_exc, :);
            
            for i_spike = 1:length(spike_indices_for_neuron)
                spike_idx = spike_indices_for_neuron(i_spike);
                x_time = all_spike_times(spike_idx);
                y_center = n_neuron;
                
                % All spikes have the same height now
                current_tick_h = tick_base_height_val;
                y_coords = [y_center - current_tick_h/2, y_center + current_tick_h/2];
                
                plot([x_time x_time], y_coords, ...
                     'Color', neuron_color_rgb, ...
                     'LineWidth', tick_line_width);
            end
        end
    end
    
    hold off;
    ylabel('Spike Raster');
    ylim([0 n+0.5]);
    yticks([])
    % set(gca, 'XColor', 'none');
    xlabel('Time (s)');
    
elseif strcmpi(sr_or_poisson, 'sr_stacked')
    % New stacked continuous spike rate plot
    hold on;
    
    % y_spacing is a factor of the 95th percentile of r_ts. This can be adjusted.
    y_spacing_factor = 3;
    p95_rate = prctile(r_ts(:), 97);
    y_spacing = y_spacing_factor * p95_rate;
    if y_spacing == 0
        y_spacing = max(r_ts(:)); % Fallback if 95th percentile is 0
        if y_spacing == 0
            y_spacing = 1; % Fallback if all rates are 0
        end
    end

    inh_color_rgb = [0.75 0 0]; % Red for inhibitory

    % Plot inhibitory neurons
    if ~isempty(I_indices)
        for k = 1:numel(I_indices)
            i = I_indices(k);
            y_level = -numel(I_indices)*y_spacing + (k-1)*y_spacing;
            plot_data = r_ts(i, plot_indices) + y_level;
            plot(t_display, plot_data, 'Color', inh_color_rgb);
        end
    end

    % Get colormap for and plot excitatory neurons
    if ~isempty(E_indices)
        cmap_exc = lines(numel(E_indices));
        % for k = 1:numel(E_indices)
        for k = numel(E_indices):-1:1
            i = E_indices(k);
            plot_color = cmap_exc(k, :);
            y_level = (i-1) * y_spacing;
            plot_data = r_ts(i, plot_indices) + y_level;
            plot(t_display, plot_data, 'Color', plot_color);
        end
    end
    
    hold off;
    ylabel('Spike Rate');
    yticks([]);
    % set(gca, 'XColor', 'none');
    xlabel('Time, s')

    % add a vertical scale bar
    scale_bar_height_data = 25;
    % Defer the call by storing the necessary info
    scale_bar_info.ax = gca;
    scale_bar_info.height = scale_bar_height_data;
    scale_bar_info.label = [num2str(scale_bar_height_data) ' Hz'];

else
    % Original spike rate plot
    hold on;
    if ~isempty(I_indices)
        plot(t_display, r_ts(I_indices, plot_indices)', 'r');
    end
    if ~isempty(E_indices)
        set(gca,'ColorOrderIndex',1); 
        plot(t_display, r_ts(E_indices, plot_indices)');
    end
    hold off;
    ylabel('Spike Rate');
end

box off;
% set(gca, 'XTickLabel', []);
% set(gca, 'XTick', []);

% % Subplot 4: SFA sum (was subplot 3)
% ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
% sp_idx = sp_idx + 1;
% a_sum_plot = zeros(n, nt); % Initialize as n x nt
% if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
%     % a_E_ts is n_E x n_a_E x nt. sum across 2nd dim -> n_E x 1 x nt. Squeeze -> n_E x nt
%     a_sum_plot(E_indices, :) = squeeze(sum(a_E_ts, 2));
% end
% if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
%     a_sum_plot(I_indices, :) = squeeze(sum(a_I_ts, 2));
% end

% hold on;
% % Plot I neurons SFA sum if they have SFA states and c_SFA is non-zero for them
% if ~isempty(I_indices) && params.n_a_I > 0
%     active_I_sfa = params.c_SFA(I_indices) ~= 0;
%     if any(active_I_sfa)
%          inh_color_rgb = [1 0 0];  % Red for inhibitory
%          plot(t_display, a_sum_plot(I_indices(active_I_sfa), plot_indices)', 'Color', inh_color_rgb);
%     end
% end
% % Plot E neurons SFA sum
% if ~isempty(E_indices) && params.n_a_E > 0
%     active_E_sfa = params.c_SFA(E_indices) ~= 0;
%     if any(active_E_sfa)
%         cmap_exc = lines(numel(E_indices));
%         for k = numel(E_indices):-1:1
%             if active_E_sfa(k)
%                 i = E_indices(k);
%                 plot_color = cmap_exc(k, :);
%                 plot(t_display, a_sum_plot(i, plot_indices)', 'Color', plot_color);
%             end
%         end
%     end
% end
% hold off;
% if params.n_a_E <= 1 && params.n_a_I <= 1
%     ylabel({'Spike Freq. Adapt., a'});
% else
%     % ylabel({'Spike Freq. Adapt.', '$\sum\limits_k a_k$'}, 'Interpreter', 'latex');
%     ylabel({'Spike Freq.', 'Adaptation'});
% end
% box off;
% set(gca, 'XTickLabel', []);
% set(gca, 'XTick', []);

% % Subplot 5: STD product (was subplot 4)
% ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
% sp_idx = sp_idx + 1;
% b_prod_plot = ones(n, nt); % Initialize as n x nt, default product is 1
% if params.n_E > 0 && params.n_b_E > 0 && ~isempty(b_E_ts)
%     % b_E_ts is n_E x n_b_E x nt. prod across 2nd dim -> n_E x 1 x nt. Squeeze -> n_E x nt
%     b_prod_plot(E_indices, :) = squeeze(prod(b_E_ts, 2));
% end
% if params.n_I > 0 && params.n_b_I > 0 && ~isempty(b_I_ts)
%     b_prod_plot(I_indices, :) = squeeze(prod(b_I_ts, 2));
% end

% hold on;
% % Plot I neurons STD product if they have STD states and F_STD is non-zero
% if ~isempty(I_indices) && params.n_b_I > 0
%     active_I_std = params.F_STD(I_indices) ~= 0;
%     if any(active_I_std)
%         inh_color_rgb = [1 0 0];  % Red for inhibitory
%         plot(t_display, b_prod_plot(I_indices(active_I_std), plot_indices)', 'Color', inh_color_rgb);
%     end
% end
% % Plot E neurons STD product
% if ~isempty(E_indices) && params.n_b_E > 0
%     active_E_std = params.F_STD(E_indices) ~= 0;
%     if any(active_E_std)
%         cmap_exc = lines(numel(E_indices));
%         for k = numel(E_indices):-1:1
%             if active_E_std(k)
%                 i = E_indices(k);
%                 plot_color = cmap_exc(k, :);
%                 plot(t_display, b_prod_plot(i, plot_indices)', 'Color', plot_color);
%             end
%         end
%     end
% end
% hold off;
% if params.n_b_E <= 1 && params.n_b_I <= 1
%     % ylabel({'Syn. Dep., b'});
%     ylabel({'Short-Term', 'Syn. Dep.'})
% else
%     ylabel({'Syn. Dep.','$\prod\limits_m b_m$'}, 'Interpreter', 'latex');
% end
% box off;
% ylim([0 1.1]); 
% yticks([0 1]);
% if strcmpi(Lya_method, 'none')
%     xlabel('Time (s)');
% else
%     set(gca, 'XTickLabel', []);
%     set(gca, 'XTick', []);
% end

% % Subplot 6: Lyapunov Exponents (was subplot 5, if calculated)
% if ~strcmpi(Lya_method, 'none')
%     ax_handles(end+1) = subplot(num_subplots, 1, sp_idx);
%     sp_idx = sp_idx + 1;
%     hold on;
%     if strcmpi(Lya_method, 'benettin')
%         % Check if there is any Lyapunov data to plot.
%         % t_lya would be empty if divergence occurred before the first interval.
%         if ~isempty(t_lya)

%             plot([t_lya(1) t_lya(end)], [LLE LLE], 'Color',[0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', sprintf('Global LLE: %.2f', LLE));
%             % plot(t_lya, finite_lya, 'r', 'LineWidth', 2, 'DisplayName', 'Finite-time LLE');
%             plot(t_lya, local_lya, 'b', 'LineWidth', 1.5, 'DisplayName', 'Local LLE');
            
%             % Adjust y-limits to ensure all data is visible
%             ylim_current = get(gca, 'YLim');
%             ylim_new = [min(ylim_current(1), LLE - 0.05*abs(LLE)), max(ylim_current(2), LLE + 0.05*abs(LLE))];
%             if all(isfinite(ylim_new))
%                 ylim(ylim_new);
%             end
%         else
%             text(0.5, 0.5, 'No Lyapunov data (diverged early)', 'Units', 'normalized', 'HorizontalAlignment', 'center');
%         end

%     elseif strcmpi(Lya_method, 'qr')
%         if ~isempty(LE_spectrum) && ~isempty(t_lya) && ~isempty(N_sys_eqs_lya)
%             colors = lines(N_sys_eqs_lya);
%             if ~isempty(local_LE_spectrum_t)
%                 for i = 1:N_sys_eqs_lya
%                     plot(t_lya, local_LE_spectrum_t(:,i), '--', 'Color', colors(i,:), 'DisplayName', sprintf('Local LE(%d)', i));
%                 end
%             end
%             if ~isempty(finite_LE_spectrum_t)
%                 for i = 1:N_sys_eqs_lya
%                     plot(t_lya, finite_LE_spectrum_t(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Finite LE(%d)', i));
%                 end
%             end
%             for i = 1:N_sys_eqs_lya
%                 plot([t_lya(1) t_lya(end)], [LE_spectrum(i) LE_spectrum(i)], ':', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('Global LE(%d): %.4f', i, LE_spectrum(i)));
%             end
%             % For QR method, fix y-limits for better comparability of spectra
%             ylim([-0.5 0.5]);
%         else
%              text(0.5, 0.5, 'No Lyapunov data to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center');
%         end
%     end
%     hold off;
%     ylabel({'Lyapunov', 'Exponent'});
%     xlabel('Time (s)');
%     % legend('show', 'Location', 'best','FontSize',10);
%     leg = legend('show', 'Location', 'North');
%     leg.Box = 'off';
%     box off;
%     ylim([-0.15 0.15])
% end

% Link all subplots' x-axes and set limits
if numel(ax_handles) > 0 && all(isgraphics(ax_handles))
    linkaxes(ax_handles, 'x');
    % Set xlim on the first axes, and linkaxes will propagate it to the others.
    axes(ax_handles(1)); 
    xlim(T_plot_limits);
end

% Add the scale bar now, after all subplots have been created
if ~isempty(scale_bar_info)
    add_scalebar_outside_subplot(scale_bar_info.ax, scale_bar_info.height, scale_bar_info.label);
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