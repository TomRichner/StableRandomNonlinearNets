% plot_parameter_space_results.m
%
% This script loads the fully consolidated data from a parameter space search
% and creates histogram plots to visualize the distribution of key metrics
% (LLE, mean spike rate) for each adaptation condition.

close all;
clear all;
clc;

%% --- Script Configuration ---

% Define whether to plot raw counts or probability density function (PDF)
% Options: 'counts' or 'pdf'
counts_or_pdf = 'pdf'; % User can change this to 'pdf'
invisible_y_axis = true;

% Define LLE histogram parameters
% lle_range = [-1.5, 1.5];
lle_range = [-1.5, 1.5];
n_bins_lle = 25;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins_lle), inf];

% Define Mean Firing Rate histogram parameters
rate_range = [0, 10];
n_bins_rate = 25;
% For rate, we use a lower bound of 0 instead of -inf
% We use n_bins_rate intervals, so we need n_bins_rate+1 edges for linspace
rate_bins = [linspace(rate_range(1), rate_range(2), n_bins_rate + 1), inf];

% Define width multiplier for outer "infinite" bins
outer_bin_width_multiplier = 1;

% Define manual ticks for the plots
% middle_ticks_lle = [-1.0, 0, 1.0];
middle_ticks_lle = [0];
middle_ticks_rate = [0, 5];

% Define Mutual Information histogram parameters
mi_delay_for_histogram_samples = 200; % Delay in samples for the MI histogram column
mi_range = [0, 0.3];
n_bins_mi = 25;
mi_bins = [linspace(mi_range(1), mi_range(2), n_bins_mi + 1), inf];
middle_ticks_mi = [0];

% Readable titles for conditions
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});
    
%% Get the base directory for the analysis
fprintf('Please select the main experiment directory to plot...\n');
base_dir = uigetdir(pwd, 'Select Experiment Directory');
if isequal(base_dir, 0)
    fprintf('User cancelled. Exiting.\n');
    return;
end
fprintf('Selected directory: %s\n', base_dir);

%% Find and load the data for each condition
conditions_to_plot = {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'};
num_conditions = length(conditions_to_plot);
all_data = struct();
found_all = true;

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    data_file = fullfile(base_dir, condition_name, sprintf('param_space_results_%s.mat', condition_name));
    
    if exist(data_file, 'file')
        fprintf('Loading: %s\n', data_file);
        all_data.(condition_name) = load(data_file);
    else
        fprintf('Error: Data file not found for condition "%s".\nExpected at: %s\n', condition_name, data_file);
        found_all = false;
    end
end

if ~found_all
    fprintf('Could not find all required data files. Exiting.\n');
    return;
end

%% Process the data to extract LLE and mean rates
extracted_values = struct();
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    results = all_data.(condition_name).results;
    
    lles = [];
    mean_rates = [];
    mi_for_hist = [];
    all_mi_vectors = {};
    shift_vec_stored = [];
    fs_stored = [];
    
    num_results = length(results);
    for j = 1:num_results
        res = results{j};
        if isstruct(res) && isfield(res, 'success') && res.success
            if isempty(fs_stored) && isfield(res, 'parameters') && isfield(res.parameters, 'fs')
                fs_stored = res.parameters.fs;
            end
            if isfield(res.simulation_output, 'LLE')
                lles(end+1) = res.simulation_output.LLE;
            end
            if isfield(res.simulation_output, 'mean_rate')
                mean_rates(end+1) = res.simulation_output.mean_rate;
            end
            
            % Process Mutual Information data
            has_mi_fields = isfield(res.simulation_output, 'MI_vs_shift') && isfield(res.simulation_output, 'shift_vec');
            
            is_valid_mi = false;
            if has_mi_fields
                shift_vec = res.simulation_output.shift_vec;
                mi_vec = res.simulation_output.MI_vs_shift;
                % Check for vector type, no NaNs, and matching length
                if isvector(mi_vec) && ~any(isnan(mi_vec)) && isvector(shift_vec) && ~any(isnan(shift_vec)) && (length(mi_vec) == length(shift_vec))
                    is_valid_mi = true;
                end
            end

            if is_valid_mi
                % This is a valid MI run
                all_mi_vectors{end+1} = mi_vec;
                
                % Find MI at the specified delay for the histogram
                [~, delay_idx] = min(abs(shift_vec - mi_delay_for_histogram_samples));
                if ~isempty(delay_idx)
                    mi_for_hist(end+1) = mi_vec(delay_idx);
                else
                    mi_for_hist(end+1) = 0; % Fallback if delay not found
                end

                % Store the template shift_vec if it's the first time
                if isempty(shift_vec_stored)
                    shift_vec_stored = shift_vec;
                end
            else
                % This run has invalid, NaN, or missing MI data. Treat as 0 MI.
                % We can only do this if we have a template for the vector size.
                if ~isempty(shift_vec_stored)
                    all_mi_vectors{end+1} = zeros(size(shift_vec_stored));
                    mi_for_hist(end+1) = 0;
                end
                % Note: If the first few runs have invalid MI, they are skipped.
                % This is acceptable as they cannot provide the template vector size.
            end
        end
    end
    extracted_values.(condition_name).lles = lles;
    extracted_values.(condition_name).mean_rates = mean_rates;
    extracted_values.(condition_name).mi_for_hist = mi_for_hist;
    extracted_values.(condition_name).all_mi_vectors = all_mi_vectors;
    extracted_values.(condition_name).shift_vec = shift_vec_stored;
    extracted_values.(condition_name).fs = fs_stored;
    fprintf('Condition "%s": Found %d LLEs, %d mean rates, %d MIs for histogram, and %d full MI vectors.\n', ...
        condition_name, length(lles), length(mean_rates), length(mi_for_hist), length(all_mi_vectors));
end

%% Create the plots
% Create a single figure with num_conditions rows x 3 columns
fig_main = figure('Name', 'LLE, Rate, and MI Distributions', 'Position', [100, 100, 810, 1020]);

% --- LLE Plots (Column 1) ---
ax_lle = gobjects(num_conditions, 1); % Store axes handles for linking
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_lle(i) = subplot(num_conditions, 3, 3*i-2); % Column 1, row i
    
    % Get histogram counts using the original infinite bins
    [counts, ~] = histcounts(extracted_values.(condition_name).lles, lle_bins);
    total_samples = sum(counts);
    
    % Create finite bins for plotting with wider outer bins
    std_width_lle = (lle_range(2) - lle_range(1)) / (n_bins_lle - 1);
    middle_edges_lle = linspace(lle_range(1), lle_range(2), n_bins_lle);
    % Correctly construct finite bins to match counts
    first_edge = middle_edges_lle(1) - outer_bin_width_multiplier * std_width_lle;
    last_edge = middle_edges_lle(end) + outer_bin_width_multiplier * std_width_lle;
    lle_bins_finite = [first_edge, middle_edges_lle, last_edge];
    
    if strcmpi(counts_or_pdf, 'pdf')
        % Calculate PDF values. The total area will be 1.
        bin_widths = diff(lle_bins_finite);
        % plot_values = counts ./ (n_bins_lle * total_samples * bin_widths);
        % plot_values = counts ./ (sum(counts) .* bin_widths);
        plot_values = counts ./ (sum(counts));

        y_label_text = 'Probability Density';
    else % Default to counts
        plot_values = counts;
        y_label_text = 'Networks';
    end

    % Plot using histogram object with specified values and finite bins
    h = histogram('BinEdges', lle_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

    % Add x-label only to the bottom row
    if i == num_conditions
        xlabel('LLE (\lambda_1)','FontSize',22);
    end
    % Add y-label with condition name to each row
    if i == 1
        % Always use text for ylabel
        text(-0.2, 0.5, sprintf('%s\n%s', condition_titles(condition_name), ''), ...
            'Units', 'normalized', 'Rotation', 90, 'FontSize', 22, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        text(-0.18, 0.5, sprintf('%s\n%s', '', y_label_text), ...
            'Units', 'normalized', 'Rotation', 90, 'FontSize', 18, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
        % Always use text for ylabel
        text(-0.2, 0.5, sprintf('%s\n%s', condition_titles(condition_name), ''), ...
            'Units', 'normalized', 'Rotation', 90, 'FontSize', 22, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    % grid on;
    box off;
    
    % Make y-axis invisible if requested
    if invisible_y_axis
        set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none');
    end
    
    % --- Custom Tick Labeling ---
    ax = gca;
    xlim([lle_bins_finite(1), lle_bins_finite(end)]);

    % Calculate center of outer bins for special labels
    center_first_bin = (lle_bins_finite(1) + lle_bins_finite(2)) / 2;
    center_last_bin = (lle_bins_finite(end-1) + lle_bins_finite(end)) / 2;
    
    final_ticks = [center_first_bin, middle_ticks_lle, center_last_bin];
    
    % Create labels for the ticks
    new_labels = cell(1, length(final_ticks));
    new_labels{1} = sprintf('<= %.1f', lle_range(1));
    for k = 1:length(middle_ticks_lle)
        new_labels{1+k} = sprintf('%.1f', middle_ticks_lle(k));
    end
    new_labels{end} = sprintf('>= %.1f', lle_range(2));
    
    xticks(ax, final_ticks);
    xticklabels(ax, new_labels);
    ax.XTickLabelRotation = 45;
end

% --- Mean Rate Plots (Column 2) ---
ax_rate = gobjects(num_conditions, 1); % Store axes handles for linking
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_rate(i) = subplot(num_conditions, 3, 3*i-1); % Column 2, row i
    
    % Get histogram counts using the original bins
    [counts, ~] = histcounts(extracted_values.(condition_name).mean_rates, rate_bins);
    total_samples = sum(counts);

    % Create finite bins for plotting with a wider outer bin
    std_width_rate = (rate_range(2) - rate_range(1)) / n_bins_rate;
    middle_edges_rate = linspace(rate_range(1), rate_range(2), n_bins_rate + 1);
    
    % Correctly construct finite bins to match counts
    first_edge = middle_edges_rate(1); % This is 0
    last_edge = middle_edges_rate(end) + outer_bin_width_multiplier * std_width_rate;
    rate_bins_finite = [first_edge, middle_edges_rate(2:end), last_edge];
    
    if strcmpi(counts_or_pdf, 'pdf')
        % Calculate PDF values. The total area will be 1.
        bin_widths = diff(rate_bins_finite);
        % plot_values = counts ./ (n_bins_rate * total_samples * bin_widths);
        % plot_values = counts ./ (sum(counts) .* bin_widths);
        plot_values = counts ./ (sum(counts));
        y_label_text = 'Probability Density';
    else % Default to counts
        plot_values = counts;
        y_label_text = 'Networks';
    end

    % Plot using histogram object with specified values and finite bins
    h = histogram('BinEdges', rate_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

    % Add x-label only to the bottom row
    if i == num_conditions
        xlabel('Mean Firing Rate (Hz)','FontSize',22);
    end
    
    % grid on;
    box off;
    
    % Make y-axis invisible if requested
    if invisible_y_axis
        set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none');
    end
    
    % --- Custom Tick Labeling ---
    ax = gca;
    xlim([rate_bins_finite(1), rate_bins_finite(end)]);
    
    % Calculate center of last bin for special label
    center_last_bin = (rate_bins_finite(end-1) + rate_bins_finite(end)) / 2;

    final_ticks = [middle_ticks_rate, center_last_bin];
    
    % Create labels for the ticks
    new_labels = cell(1, length(final_ticks));
    for k = 1:length(middle_ticks_rate)
        new_labels{k} = sprintf('%.0f', middle_ticks_rate(k));
    end
    new_labels{end} = sprintf('>= %.0f', rate_range(2));

    xticks(ax, final_ticks);
    xticklabels(ax, new_labels);
    ax.XTickLabelRotation = 45;
end

% --- MI Plots (Column 3) ---
ax_mi = gobjects(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_mi(i) = subplot(num_conditions, 3, 3*i); % Column 3, row i
    
    % Get histogram counts
    [counts, ~] = histcounts(extracted_values.(condition_name).mi_for_hist, mi_bins);
    total_samples = sum(counts);

    % Create finite bins for plotting
    std_width_mi = (mi_range(2) - mi_range(1)) / n_bins_mi;
    middle_edges_mi = linspace(mi_range(1), mi_range(2), n_bins_mi + 1);
    first_edge = middle_edges_mi(1);
    last_edge = middle_edges_mi(end) + outer_bin_width_multiplier * std_width_mi;
    mi_bins_finite = [first_edge, middle_edges_mi(2:end), last_edge];
    
    if strcmpi(counts_or_pdf, 'pdf')
        bin_widths = diff(mi_bins_finite);
        plot_values = counts ./ (sum(counts));
        y_label_text = 'Probability Density';
    else % Default to counts
        plot_values = counts;
        y_label_text = 'Networks';
    end

    h = histogram('BinEdges', mi_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

    if i == num_conditions
        actual_delay_for_label = mi_delay_for_histogram_samples;
        current_shift_vec = extracted_values.(condition_name).shift_vec;
        if ~isempty(current_shift_vec)
            [~, idx] = min(abs(current_shift_vec - mi_delay_for_histogram_samples));
            actual_delay_for_label = current_shift_vec(idx);
        end
        xlabel(sprintf('MI (bits)', actual_delay_for_label),'FontSize',22);
    end
    
    box off;
    
    if invisible_y_axis
        set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none');
    end
    
    % --- Custom Tick Labeling ---
    ax = gca;
    xlim([mi_bins_finite(1), mi_bins_finite(end)]);
    
    center_last_bin = (mi_bins_finite(end-1) + mi_bins_finite(end)) / 2;

    final_ticks = [middle_ticks_mi, center_last_bin];
    
    new_labels = cell(1, length(final_ticks));
    for k = 1:length(middle_ticks_mi)
        new_labels{k} = sprintf('%.1f', middle_ticks_mi(k));
    end
    new_labels{end} = sprintf('>= %.1f', mi_range(2));

    xticks(ax, final_ticks);
    xticklabels(ax, new_labels);
    ax.XTickLabelRotation = 45;
end

% Link all y-axes together
all_axes = [ax_lle; ax_rate; ax_mi];
linkaxes(all_axes, 'y');
% linkaxes(ax_lle, 'y');
% linkaxes(ax_rate, 'y');

% %% Add letters to subplots
% num_subplots = 2 * num_conditions;
% if num_subplots > 0
%     if num_subplots <= 26
%         letters = arrayfun(@(x) sprintf('(%c)', x), 'a':char('a' + num_subplots - 1), 'UniformOutput', false);
%         AddLetters2Plots(fig_main, letters, 'FontSize', 20, 'FontWeight', 'Normal', 'HShift', -0.027, 'VShift', -0.03, 'Location', 'NorthWest');
%     else
%         error('More than 26 subplots, cannot add panel letters.');
%     end
% end

%% Add a letter to the first subplot of each row
num_subplots = num_conditions * 3;

letters = cell(1, num_subplots);
letter_char_code = 'a';
current_subplot_idx = 1;

% Iterate through rows columns to build
% the cell array of labels.
for i_row = 1:num_conditions
    for i_col = 1:3
        if i_col == 1 % Only add a letter to the first column
            letters{current_subplot_idx} = sprintf('(%c)', letter_char_code);
            letter_char_code = letter_char_code + 1;
        else
            letters{current_subplot_idx} = '';
        end
        current_subplot_idx = current_subplot_idx + 1;
    end
end

AddLetters2Plots(fig_main, letters, 'FontSize', 20, 'FontWeight', 'Normal', 'HShift', -0.026, 'VShift', -0.028, 'Location', 'NorthWest');

%% --- Create a separate figure for MI Dynamics ---
fig_mi_dyn = figure('Name', 'MI Dynamics vs Delay', 'Position', [150, 150, 400, 900]);
ax_mi_dyn = gobjects(num_conditions, 1);
global_max_count = 0; % To store the absolute max count across all subplots

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_mi_dyn(i) = subplot(num_conditions, 1, i);
    
    all_mi_vectors = extracted_values.(condition_name).all_mi_vectors;
    shift_vec = extracted_values.(condition_name).shift_vec;
    fs = extracted_values.(condition_name).fs;

    if isempty(all_mi_vectors) || isempty(shift_vec) || isempty(fs)
        title(condition_titles(condition_name));
        text(0.5, 0.5, 'No MI Data', 'HorizontalAlignment', 'center');
        axis off;
        continue;
    end

    shift_vec_sec = shift_vec / fs;

    % Create the 2D histogram matrix slice by slice
    mi_matrix = vertcat(all_mi_vectors{:});
    y_edges = linspace(mi_range(1), mi_range(2), n_bins_mi + 1);
    mi_hist_matrix = zeros(n_bins_mi, length(shift_vec));

    for i_delay = 1:length(shift_vec)
        mi_values_at_this_delay = mi_matrix(:, i_delay);
        counts = histcounts(mi_values_at_this_delay, y_edges);
        mi_hist_matrix(:, i_delay) = counts(:);
    end
    
    % Update the global maximum count
    current_max = max(mi_hist_matrix, [], 'all');
    if current_max > global_max_count
        global_max_count = current_max;
    end

    % Plot the 2D histogram
    y_bin_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    imagesc(shift_vec_sec, y_bin_centers, mi_hist_matrix);
    axis xy;
    colormap(hot);
    
    title(condition_titles(condition_name), 'FontSize', 22);
    
    if i == num_conditions
        xlabel('Delay (s)','FontSize',22);
    end
    ylabel('MI (bits)', 'FontSize', 22);

    % Make sure ticks are reasonable for x-axis
    xticks(linspace(shift_vec_sec(1), shift_vec_sec(end), 5));
    xtickformat('%.2f');
end

% Apply the unified color limit to all subplots
for i = 1:num_conditions
    if ishandle(ax_mi_dyn(i)) && ~isempty(ax_mi_dyn(i).Children)
         caxis(ax_mi_dyn(i), [0, 0.9*global_max_count]);
    end
end

linkaxes(ax_mi_dyn, 'y');

%% --- Create a third figure for Median MI vs. Delay ---
fig_mi_median = figure('Name', 'Median MI vs. Delay', 'Position', [200, 200, 400, 900]);
ax_mi_median = gobjects(num_conditions, 1);

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_mi_median(i) = subplot(num_conditions, 1, i);
    
    all_mi_vectors = extracted_values.(condition_name).all_mi_vectors;
    shift_vec = extracted_values.(condition_name).shift_vec;
    fs = extracted_values.(condition_name).fs;

    if isempty(all_mi_vectors) || isempty(shift_vec) || isempty(fs)
        title(condition_titles(condition_name));
        text(0.5, 0.5, 'No MI Data', 'HorizontalAlignment', 'center');
        axis off;
        continue;
    end
    
    % Convert MI vectors cell array to a matrix
    mi_matrix = vertcat(all_mi_vectors{:});
    
    % Calculate the median MI at each delay
    median_mi_vs_delay = median(mi_matrix, 1, 'omitnan');
    
    % Convert delay vector to seconds
    shift_vec_sec = shift_vec / fs;
    
    % Plot the median line
    plot(shift_vec_sec, median_mi_vs_delay, 'k-', 'LineWidth', 2);
    
    title(condition_titles(condition_name), 'FontSize', 22);
    box off;
    grid on;
    
    if i == num_conditions
        xlabel('Delay (s)','FontSize',22);
    end
    ylabel('Median MI (bits)', 'FontSize', 22);

    % Make sure ticks are reasonable for x-axis
    xticks(linspace(shift_vec_sec(1), shift_vec_sec(end), 5));
    xtickformat('%.2f');
end

linkaxes(ax_mi_median, 'y');


%% Save figures
output_dir_for_figs = fullfile(base_dir, 'analysis_plots');
if ~exist(output_dir_for_figs, 'dir')
    mkdir(output_dir_for_figs);
end

fprintf('\nSaving figures to: %s\n', output_dir_for_figs);
save_some_figs_to_folder_2(output_dir_for_figs, 'LLE_rate_and_MI_distributions_v1', fig_main.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'MI_dynamics_vs_delay_v1', fig_mi_dyn.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'Median_MI_vs_delay_v1', fig_mi_median.Number, {'fig', 'svg', 'png'});

fprintf('Plotting complete.\n');
beep; pause(0.5); beep % wake up