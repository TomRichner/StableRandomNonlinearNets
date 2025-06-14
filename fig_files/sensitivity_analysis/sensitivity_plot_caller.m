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
custom_x_labels('n') = '# Neurons in Network';

% Define scaling factors for x-tick labels for specific parameters
rescale_x_ticks = containers.Map('KeyType', 'char', 'ValueType', 'double');
rescale_x_ticks('EE_factor') = 0.5; % Example: if EE_factor is varied, scale ticks by 0.5
% Add other parameters here if they need rescaling, e.g.,
% rescale_x_ticks('some_other_param') = 10;

% --- NEW: readable titles for each condition -----------------------------
custom_condition_titles = containers.Map('KeyType','char','ValueType','char');
custom_condition_titles('no_adaptation') = 'No Adaptation';
custom_condition_titles('sfa_only')      = 'SFA Only';
custom_condition_titles('std_only')      = 'STD Only';
custom_condition_titles('sfa_and_std')   = 'SFA + STD';
% -------------------------------------------------------------------------

%% Set up paths and parameters
% -------------------------------------------------------------------------
% Always work relative to the folder that contains THIS script so that the
% user's current working directory does not matter.
% -------------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));   % Absolute path of this .m file

% sensitivity_results_dir = 'sensi_c_SFA_and_tau_b_E_reduced_nLevs_9_nReps_16_jun_12_25_ 3_31_pm';
% sensitivity_results_dir = 'sensi_longer_tau_a_E_nLevs_9_nReps_16_jun_12_25_ 4_28_pm';
% sensitivity_results_dir = 'sensi_n_b_E_is_2_nLevs_9_nReps_16_jun_12_25_ 4_37_pm';
% sensitivity_results_dir = 'sensi_multi_nLevs_5_nReps_6_jun_12_25_ 4_56_pm';
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_5_nReps_6_jun_12_25_ 5_09_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_9_nReps_16_jun_12_25_ 5_19_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_21_nReps_100_jun_12_25_10_53_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_11_nReps_25_jun_13_25_ 8_15_am'
% sensitivity_results_dir = 'sensi_n_b_2_tau_a_E_2_15_c_SFA_0p5_nLevs_11_nReps_25_jun_13_25_ 8_47_am'
% sensitivity_results_dir = 'sensi_n_b_2_tau_a_E_2_15_c_SFA_0p5_nLevs_15_nReps_25_jun_13_25_ 7_35_pm'
% sensitivity_results_dir = 'sensi_n_b_2_tau_a_E_2_15_c_SFA_0p5_nLevs_15_nReps_25_jun_13_25_ 9_50_pm';
% sensitivity_results_dir = 'sensi_quick_LLE_and_SR_nLevs_5_nReps_5_jun_13_25_11_57_pm'
% sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_n_nLevs_25_nReps_25_jun_14_25_12_25_am'
% sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_nLevs_7_nReps_10_jun_14_25_ 6_53_am'
% sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_nLevs_7_nReps_10_jun_14_25_ 7_05_am';
% sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_nLevs_13_nReps_25_jun_14_25_ 7_12_am';
% sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_randWalk_nLevs_25_nReps_50_jun_14_25_ 4_55_pm'
sensitivity_results_dir = 'sensi_LLE_and_SR_EE_W_randWalk_nLevs_25_nReps_50_jun_14_25_10_04_pm';

% If the user left the path relative, prepend the script directory
if ~isfolder(sensitivity_results_dir)
    sensitivity_results_dir = fullfile(script_dir, sensitivity_results_dir);
end

output_dir_base = sensitivity_results_dir;

% Create base output directory if it doesn't exist
if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

% Define LLE histogram parameters
lle_range = [-1.5, 1.5];
n_bins = 25;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins), inf];

% Define Mean Firing Rate histogram parameters
rate_range = [0, 5];
n_bins_rate = 25;
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

% Create a separate large figure for each metric
main_fig_lle = figure('Name', 'LLE Sensitivity', 'Position', [100, 100, 250 * n_conditions, 250 * n_params], 'Visible', 'on');
main_fig_rate = figure('Name', 'Mean Rate Sensitivity', 'Position', [150, 150, 250 * n_conditions, 250 * n_params], 'Visible', 'on');

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
        subplot_idx = (i - 1) * n_conditions + c_idx;
        sp_ax_lle = subplot(n_params, n_conditions, subplot_idx);
        plot_single_metric(sp_ax_lle, param_file, lle_bins, 'LLE', '$\lambda_1$', ...
                           param_name, condition_name, condition_title, ...
                           x_label_str, scale_factor, i, c_idx);
        
        % --- Mean Rate Plot ---
        figure(main_fig_rate);
        sp_ax_rate = subplot(n_params, n_conditions, subplot_idx);
        plot_single_metric(sp_ax_rate, param_file, rate_bins, 'mean_rate', 'Mean Rate (Hz)', ...
                           param_name, condition_name, condition_title, ...
                           x_label_str, scale_factor, i, c_idx);
    end
end

% Save the combined figures
save_some_figs_to_folder_2(output_dir_base, 'sensitivity_LLE_comparison_all_params', main_fig_lle.Number, {'png', 'svg', 'fig'});
save_some_figs_to_folder_2(output_dir_base, 'sensitivity_rate_comparison_all_params', main_fig_rate.Number, {'png', 'svg', 'fig'});

% close(main_fig_lle); % Don't close figure after saving to allow inspection
% close(main_fig_rate);

fprintf('\n=== Sensitivity plotting complete ===\n');
fprintf('Plots saved to: %s\n', output_dir_base); 


function plot_single_metric(sp_ax, param_file, hist_bins, metric_name, y_label_metric, ...
                            param_name, condition_name, condition_title, ...
                            x_label_str, scale_factor, i, c_idx)
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
    
    try
        % Call the plotting function which returns a handle to an invisible figure
        temp_fig_handle = sensitivity_plot(param_file, hist_bins, metric_name, y_label_metric);
        
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

        % Manually copy colorbar
        temp_cb = findobj(temp_fig_handle, 'type', 'colorbar');
        if ~isempty(temp_cb)
            cb_new = colorbar(sp_ax);
            set(cb_new, 'Limits', get(temp_cb, 'Limits'));
        end
        
        % Set titles and labels for the grid
        if i == 1           % top-row → write condition title
            title(sp_ax, condition_title, ...
                  'Interpreter','none', ...   % plain text
                  'FontWeight','normal');     % not bold
        end
        
        % Set y-label for the first column only
        if c_idx == 1 
            ylabel(sp_ax, y_label_metric, 'Interpreter', 'latex', 'FontSize', 22);
        end

        % Set custom x-label for the row, assuming LaTeX interpreter
        xlabel(sp_ax, x_label_str, ...
               'Interpreter','none', ...      % plain text
               'FontWeight','normal');
        
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

    catch ME
        fprintf('Error plotting %s for condition %s: %s\n', param_name, condition_name, ME.message);
        if ~isempty(ME.stack)
            fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        % Also set labels for error plots for grid consistency
        if i == 1, title(sp_ax, strrep(condition_name, '_', ' ')); end
        if c_idx == 1, ylabel(sp_ax, strrep(param_name, '_', '\_')); end
    end
end 