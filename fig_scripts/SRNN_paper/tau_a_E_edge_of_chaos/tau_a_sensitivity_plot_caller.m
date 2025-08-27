% Script for plotting sensitivity analysis results
% This script loads the sensitivity analysis results and creates visualizations
% of the largest Lyapunov exponent (LLE) across parameter ranges
 
clear;
clc;
close all;


%% Customization for plots
% Define custom x-axis labels for specific parameters. Uses LaTeX formatting.
custom_x_labels = containers.Map('KeyType', 'char', 'ValueType', 'char');
custom_x_labels('EE_factor') = 'Mean E-to-E Weight';
custom_x_labels('mean_weight') = 'Mean Weight';
custom_x_labels('n') = '\# Neurons in Network';
custom_x_labels('tau_a_E_2') = '$\max(\tau_a),  s$';
custom_x_labels('tau_b_E_2') = '$\tau_b,  s$';

% Define scaling factors for x-tick labels for specific parameters
rescale_x_ticks = containers.Map('KeyType', 'char', 'ValueType', 'double');
rescale_x_ticks('EE_factor') = 0.5; % Example: if EE_factor is varied, scale ticks by 0.5
% Add other parameters here if they need rescaling, e.g.,
% rescale_x_ticks('some_other_param') = 10;

% Define custom y-axis ticks for specific parameters and metrics
custom_y_ticks_lle = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Example: custom_y_ticks_lle('n') = [-0.5, 0, 0.5];
custom_y_ticks_lle('n') = [-0.5, 0, 0.5];
custom_y_ticks_lle('EE_factor') = [-0.5, 0, 0.5];
custom_y_ticks_lle('mean_weight') = [-0.5, 0, 0.5];
custom_y_ticks_lle('tau_a_E_2') = [-0.3, -0.15, 0];
custom_y_ticks_lle('tau_b_E_2') = [-0.3, -0.15, 0];

custom_y_ticks_rate = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Example: custom_y_ticks_rate('n') = [0, 10, 20];x

% --- Define custom x-tick vectors for specific parameters ---
custom_x_ticks = containers.Map('KeyType', 'char', 'ValueType', 'any');
custom_x_ticks('tau_a_E_2') = [5 20 40 60];
custom_x_ticks('tau_b_E_2') = [5 20 40 60];

% --- Readable titles for each condition -----------------------------
custom_condition_titles = containers.Map('KeyType','char','ValueType','char');
custom_condition_titles('no_adaptation') = 'No Adaptation';
custom_condition_titles('sfa_only')      = 'SFA Only';
custom_condition_titles('std_only')      = 'STD Only';
custom_condition_titles('sfa_and_std')   = 'SFA + STD';
% -------------------------------------------------------------------------

%% Set up paths and parameters
script_dir = fileparts(mfilename('fullpath'));   % Absolute path of this .m file

% --- Get the base directory for the analysis ---
fprintf('Please select the main sensitivity experiment directory to plot...\n');
sensitivity_results_dir = uigetdir(pwd, 'Select Experiment Directory');
if isequal(sensitivity_results_dir, 0)
    fprintf('User cancelled. Exiting.\n');
    return;
end
fprintf('Selected directory: %s\n', sensitivity_results_dir);

output_dir_base = sensitivity_results_dir;

% Create base output directory if it doesn't exist
if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

% Define LLE histogram parameters
lle_range = [-0.3, 0.1];
n_bins = 35;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins), inf];

% Define Mean Firing Rate histogram parameters
rate_range = [0, 5];
n_bins_rate = 35;
rate_bins = [linspace(rate_range(1), rate_range(2), n_bins_rate), inf];

%% Load summary file to get list of analyzed parameters and conditions
summary_file = fullfile(sensitivity_results_dir, 'sensitivity_analysis_summary_all_conditions.mat');
if ~exist(summary_file, 'file')
    error('Summary file not found: %s. Please run sensitivity_analysis.m first.', summary_file);
end

fprintf('Loading summary file: %s\n', summary_file);
summary_data_loaded = load(summary_file);
summary_data = summary_data_loaded.summary_data;

% Get conditions from the summary
if isfield(summary_data, 'conditions')
    conditions = summary_data.conditions;
    fprintf('Found %d conditions to process.\n', length(conditions));
else
    error('No conditions found in summary file. The data format may be old or corrupted.');
end

% Collect all unique parameter names from all conditions
all_param_names = {};
for c_idx = 1:length(conditions)
    condition_name = conditions{c_idx}.name;
    if isfield(summary_data, 'all_conditions_summary') && isfield(summary_data.all_conditions_summary, condition_name)
        param_names = fieldnames(summary_data.all_conditions_summary.(condition_name));
        all_param_names = [all_param_names; param_names];
    else
        fprintf('Warning: Could not find summary for condition "%s" in summary file.\n', condition_name);
    end
end
unique_param_names = unique(all_param_names);

fprintf('Found %d unique parameters to plot: %s\n', length(unique_param_names), strjoin(unique_param_names, ', '));

n_params = length(unique_param_names);
n_conditions = length(conditions);

%% Extract n_reps from the data
% Try to get n_reps from summary data first
n_reps = [];
if isfield(summary_data, 'n_reps')
    n_reps = summary_data.n_reps;
    fprintf('Found n_reps in summary data: %d\n', n_reps);
elseif isfield(summary_data, 'analysis_parameters') && isfield(summary_data.analysis_parameters, 'n_reps')
    n_reps = summary_data.analysis_parameters.n_reps;
    fprintf('Found n_reps in analysis parameters: %d\n', n_reps);
end

% If n_reps not found in summary, load it from the first available parameter file
if isempty(n_reps)
    fprintf('n_reps not found in summary data, loading from parameter file...\n');
    for c_idx = 1:n_conditions
        condition_name = conditions{c_idx}.name;
        condition_results_dir = fullfile(sensitivity_results_dir, condition_name);
        for p_idx = 1:n_params
            param_name = strtrim(unique_param_names{p_idx});
            param_file = fullfile(condition_results_dir, sprintf('sensitivity_%s.mat', param_name));
            if exist(param_file, 'file')
                fprintf('Loading n_reps from: %s\n', param_file);
                temp_data = load(param_file);
                if isfield(temp_data, 'metadata') && isfield(temp_data.metadata, 'n_reps')
                    n_reps = temp_data.metadata.n_reps;
                    fprintf('Found n_reps: %d\n', n_reps);
                    break;
                end
            end
        end
        if ~isempty(n_reps), break; end
    end
end

% Default fallback if still not found
if isempty(n_reps)
    n_reps = 25; % Default value
    fprintf('Warning: Could not find n_reps in data, using default value: %d\n', n_reps);
end

% Create a separate large figure for each metric
main_fig_lle = figure('Name', 'LLE Sensitivity', 'Position', [100, 100, 400 * n_params, 400 * n_conditions], 'Visible', 'on');
main_fig_rate = figure('Name', 'Mean Rate Sensitivity', 'Position', [150, 150, 400 * n_params, 400 * n_conditions], 'Visible', 'on');

%% Process each parameter across all conditions
for i = 1:n_params
    %––– Remove any leading / trailing blanks coming from MAT-files –––%
    param_name_raw = unique_param_names{i};
    param_name     = strtrim(param_name_raw);

    fprintf('\n\n--- Processing Parameter Row: %s ---\n', upper(param_name));

    % Get custom label and scale factor for this parameter row
    if isKey(custom_x_labels, param_name)
        x_label_str = custom_x_labels(param_name);
    else
        x_label_str = strrep(param_name, '_', '\_'); % Default label
    end

    if isKey(rescale_x_ticks, param_name)
        scale_factor = rescale_x_ticks(param_name);
    else
        scale_factor = 1.0; % Default: no scaling
    end

    % Get custom x-ticks for this parameter row
    if isKey(custom_x_ticks, param_name)
        x_ticks_vec = custom_x_ticks(param_name);
    else
        x_ticks_vec = []; % Default: auto-ticks
    end

    % Get custom y-ticks for this parameter row
    if isKey(custom_y_ticks_lle, param_name)
        y_ticks_lle = custom_y_ticks_lle(param_name);
    else
        y_ticks_lle = []; % Default: use auto-ticks
    end
    if isKey(custom_y_ticks_rate, param_name)
        y_ticks_rate = custom_y_ticks_rate(param_name);
    else
        y_ticks_rate = []; % Default: use auto-ticks
    end

    fprintf('--> Using xlabel: "%s" and x-tick scale factor: %.2f\n', x_label_str, scale_factor);

    for c_idx = 1:n_conditions
        current_condition = conditions{c_idx};
        condition_name    = current_condition.name;

        % Pick a pretty title for this condition
        if isKey(custom_condition_titles, condition_name)
            condition_title = custom_condition_titles(condition_name);
        else
            condition_title = strrep(condition_name, '_', ' ');
        end
        fprintf('--- Plotting for condition: %s → "%s" (%d/%d) ---\n', ...
                condition_name, condition_title, c_idx, n_conditions);

        condition_results_dir = fullfile(sensitivity_results_dir, condition_name);
        param_file = fullfile(condition_results_dir, ...
                             sprintf('sensitivity_%s.mat', param_name));
        
        % Verbose check so we know exactly what file we are looking for
        fprintf('      looking for: %s\n', param_file);

        % --- LLE Plot ---
        figure(main_fig_lle);
        subplot_idx = (c_idx - 1) * n_params + i;
        sp_ax_lle = subplot(n_conditions, n_params, subplot_idx);
        plot_single_metric(sp_ax_lle, param_file, lle_bins, 'LLE', '$\lambda_1$', ...
                           param_name, condition_name, condition_title, ...
                           x_label_str, scale_factor, i, c_idx, y_ticks_lle, x_ticks_vec);
        
        % --- Mean Rate Plot ---
        figure(main_fig_rate);
        sp_ax_rate = subplot(n_conditions, n_params, subplot_idx);
        plot_single_metric(sp_ax_rate, param_file, rate_bins, 'mean_rate', 'Mean Rate (Hz)', ...
                           param_name, condition_name, condition_title, ...
                           x_label_str, scale_factor, i, c_idx, y_ticks_rate, x_ticks_vec);
    end
end

%% Create separate colorbar figures
% Create colorbar figure for LLE
colorbar_fig_lle = figure('Name', 'LLE Colorbar', 'Position', [200, 200, 200, 400], 'Visible', 'on');
ax_cb_lle = axes('Position', [0.2, 0.1, 0.6, 0.8]);
cb_lle = colorbar(ax_cb_lle, 'Location', 'west');
colormap(ax_cb_lle, flipud(gray));
caxis(ax_cb_lle, [0, n_reps]);
% Set custom ticks and labels to show percentages
cb_lle.Ticks = linspace(0, n_reps, 6);  % 6 tick marks from 0 to n_reps
cb_lle.TickLabels = arrayfun(@(x) sprintf('%.1f', x/n_reps), cb_lle.Ticks, 'UniformOutput', false);
cb_lle.Label.String = 'Probability';
cb_lle.Label.FontSize = 20;
axis(ax_cb_lle, 'off');

% Create colorbar figure for Mean Rate  
colorbar_fig_rate = figure('Name', 'Mean Rate Colorbar', 'Position', [250, 250, 200, 400], 'Visible', 'on');
ax_cb_rate = axes('Position', [0.2, 0.1, 0.6, 0.8]);
cb_rate = colorbar(ax_cb_rate, 'Location', 'west');
colormap(ax_cb_rate, flipud(gray));
caxis(ax_cb_rate, [0, n_reps]);
% Set custom ticks and labels to show percentages
cb_rate.Ticks = linspace(0, n_reps, 6);  % 6 tick marks from 0 to n_reps
cb_rate.TickLabels = arrayfun(@(x) sprintf('%.1f', x/n_reps), cb_rate.Ticks, 'UniformOutput', false);
cb_rate.Label.String = 'Probability';
cb_rate.Label.FontSize = 20;
axis(ax_cb_rate, 'off');

%% add letter to each figure
num_subplots = n_params *  n_conditions;
if num_subplots > 0
    % Generate letters (A), (B), ... up to (Z)
    if num_subplots <= 26
        letters = arrayfun(@(x) sprintf('%c', x), 'A':char('A'+num_subplots-1), 'UniformOutput', false);
        AddLetters2Plots(main_fig_lle, letters, 'FontSize', 30, 'FontWeight', 'Normal', 'HShift', -0.05, 'VShift', -0.05, 'Location', 'NorthWest');
        AddLetters2Plots(main_fig_rate, letters, 'FontSize', 30, 'FontWeight', 'Normal', 'HShift', -0.05, 'VShift', -0.05, 'Location', 'NorthWest');
    else
        error('more than 26 subplots, out of letters')
    end
end

% Save the combined figures
save_some_figs_to_folder_2([output_dir_base filesep 'figs'], 'sensitivity_LLE_comparison_all_params', main_fig_lle.Number, {'png', 'svg', 'fig'});
save_some_figs_to_folder_2([output_dir_base filesep 'figs'], 'sensitivity_rate_comparison_all_params', main_fig_rate.Number, {'png', 'svg', 'fig'});

% Save the colorbar figures
save_some_figs_to_folder_2([output_dir_base filesep 'figs'], 'sensitivity_LLE_colorbar', colorbar_fig_lle.Number, {'png', 'svg', 'fig'});
save_some_figs_to_folder_2([output_dir_base filesep 'figs'], 'sensitivity_rate_colorbar', colorbar_fig_rate.Number, {'png', 'svg', 'fig'});

% close(main_fig_lle); % Don't close figure after saving to allow inspection
% close(main_fig_rate);

fprintf('\n=== Sensitivity plotting complete ===\n');
fprintf('Plots saved to: %s\n', output_dir_base); 


function plot_single_metric(sp_ax, param_file, hist_bins, metric_name, y_label_metric, ...
                            param_name, condition_name, condition_title, ...
                            x_label_str, scale_factor, i, c_idx, custom_y_ticks, custom_x_ticks)
    if ~exist(param_file, 'file')
        fprintf('Warning: File not found: %s. Skipping plot for this condition.\n', param_file);
        if i == 1
            title(strrep(condition_name, '_', ' '));
        end
        if c_idx == 1
            ylabel(sp_ax, strrep(param_name, '_', '\_'));
        end
        axis off;
        text(0.5, 0.5, 'Data not found', 'HorizontalAlignment', 'center', 'Parent', sp_ax);
        return;
    end
    
    % Call the plotting function which returns a handle to an invisible figure
    temp_fig_handle = tau_a_sensitivity_plot(param_file, hist_bins, metric_name, y_label_metric, custom_y_ticks, custom_x_ticks);
    
    % Get axes from the temporary figure
    temp_ax = findobj(temp_fig_handle, 'type', 'axes');
    
    % Copy the contents of the old axes to the new subplot axes
    copyobj(allchild(temp_ax), sp_ax);
    
    % Copy essential properties
    set(sp_ax, 'XLim', get(temp_ax, 'XLim'), 'YLim', get(temp_ax, 'YLim'), ...
        'YDir', get(temp_ax, 'YDir'), 'Colormap', get(temp_ax, 'Colormap'), ...
        'CLim', get(temp_ax, 'CLim'));

    % Copy y-tick labels from the temporary figure
    set(sp_ax, 'YTick', get(temp_ax, 'YTick'), 'YTickLabel', get(temp_ax, 'YTickLabel'));

    % Copy x-tick labels from the temporary figure to ensure they are preserved
    set(sp_ax, 'XTick', get(temp_ax, 'XTick'), 'XTickLabel', get(temp_ax, 'XTickLabel'));

    % REMOVED: Colorbar creation code
    % temp_cb = findobj(temp_fig_handle, 'type', 'colorbar');
    % if ~isempty(temp_cb)
    %     cb_new = colorbar(sp_ax);
    %     set(cb_new, 'Limits', get(temp_cb, 'Limits'));
    % end
    
    % Set y-label for the first column only
    if c_idx == 1 
        ylabel(sp_ax, y_label_metric, 'Interpreter', 'latex', 'FontSize', 22);
    end

    % Set custom x-label for the row, assuming LaTeX interpreter
    xlabel(sp_ax, x_label_str, 'Interpreter', 'latex', 'FontWeight','normal');
    
    % Rescale x-ticks if a factor is specified
    if scale_factor ~= 1.0
        drawnow; % Ensure ticks are updated before getting them
        current_ticks = xticks(sp_ax);
        new_tick_labels = arrayfun(@(x) sprintf('%.2g', x * scale_factor), current_ticks, 'UniformOutput', false);
        xticklabels(sp_ax, new_tick_labels);
    end
    
    % Close the invisible temporary figure
    close(temp_fig_handle);
    fprintf('Successfully plotted %s for condition %s\n', param_name, condition_name);
end 