% Script for plotting sensitivity analysis results
% This script loads the sensitivity analysis results and creates visualizations
% of the largest Lyapunov exponent (LLE) across parameter ranges
 
clear;
clc;
close all;


%% Set up paths and parameters
% sensitivity_results_dir = 'sensi_c_SFA_and_tau_b_E_reduced_nLevs_9_nReps_16_jun_12_25_ 3_31_pm';
% sensitivity_results_dir = 'sensi_longer_tau_a_E_nLevs_9_nReps_16_jun_12_25_ 4_28_pm';
% sensitivity_results_dir = 'sensi_n_b_E_is_2_nLevs_9_nReps_16_jun_12_25_ 4_37_pm';
% sensitivity_results_dir = 'sensi_multi_nLevs_5_nReps_6_jun_12_25_ 4_56_pm';
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_5_nReps_6_jun_12_25_ 5_09_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_9_nReps_16_jun_12_25_ 5_19_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_21_nReps_100_jun_12_25_10_53_pm'
% sensitivity_results_dir = 'sensi_multi_lower_tau_a_tau_b_nLevs_11_nReps_25_jun_13_25_ 8_15_am'
sensitivity_results_dir = 'sensi_n_b_2_tau_a_E_2_15_c_SFA_0p5_nLevs_11_nReps_25_jun_13_25_ 8_47_am'

output_dir_base = fullfile(fileparts(mfilename('fullpath')), 'plots');

% Create base output directory if it doesn't exist
if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

% Define LLE histogram parameters
lle_range = [-1.2, 1.2];
n_bins = 30;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins), inf];

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

% Create a single large figure for all plots
main_fig = figure('Position', [100, 100, 350 * n_conditions, 300 * n_params], 'Visible', 'on');
sgtitle(main_fig, 'Sensitivity Analysis Comparison', 'FontSize', 16, 'FontWeight', 'bold');

%% Process each parameter across all conditions
for i = 1:n_params
    param_name = unique_param_names{i};
    fprintf('\n\n--- Processing Parameter: %s ---\n', upper(param_name));

    for c_idx = 1:n_conditions
        current_condition = conditions{c_idx};
        condition_name = current_condition.name;
        fprintf('--- Plotting for condition: %s (%d/%d) ---\n', condition_name, c_idx, n_conditions);

        condition_results_dir = fullfile(sensitivity_results_dir, condition_name);
        param_file = fullfile(condition_results_dir, sprintf('sensitivity_%s.mat', param_name));
        
        % Calculate subplot index
        subplot_idx = (i - 1) * n_conditions + c_idx;
        sp_ax = subplot(n_params, n_conditions, subplot_idx);

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
            continue;
        end
        
        try
            % Call the plotting function which returns a handle to an invisible figure
            temp_fig_handle = sensitivity_plot(param_file, lle_bins);
            
            % Get axes from the temporary figure
            temp_ax = findobj(temp_fig_handle, 'type', 'axes');
            
            % Copy the contents of the old axes to the new subplot axes
            copyobj(allchild(temp_ax), sp_ax);
            
            % Copy essential properties
            set(sp_ax, 'XLim', get(temp_ax, 'XLim'), 'YLim', get(temp_ax, 'YLim'), ...
                'YDir', get(temp_ax, 'YDir'), 'Colormap', get(temp_ax, 'Colormap'), ...
                'CLim', get(temp_ax, 'CLim'));

            % Manually copy colorbar
            temp_cb = findobj(temp_fig_handle, 'type', 'colorbar');
            if ~isempty(temp_cb)
                cb_new = colorbar(sp_ax);
                set(cb_new, 'Limits', get(temp_cb, 'Limits'));
            end
            
            % Set titles and labels for the grid
            if i == 1 % Top row: set condition names as titles
                title(sp_ax, strrep(condition_name, '_', ' '));
            end
            
            % Set y-label for the first column only
            if c_idx == 1 
                ylabel(sp_ax, '$\lambda_1$', 'Interpreter', 'latex', 'FontSize', 14);
            end

            % Set x-label for all plots from the source plot
            xlabel(sp_ax, get(get(temp_ax, 'XLabel'), 'String'));
            
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
end

% Save the combined figure
fig_filename_png = fullfile(output_dir_base, 'sensitivity_comparison_all_params.png');
saveas(main_fig, fig_filename_png);
fig_filename_fig = fullfile(output_dir_base, 'sensitivity_comparison_all_params.fig');
saveas(main_fig, fig_filename_fig);
fprintf('Saved combined plot: %s\n', fig_filename_png);

% close(main_fig); % Don't close figure after saving to allow inspection

fprintf('\n=== Sensitivity plotting complete ===\n');
fprintf('Plots saved to: %s\n', output_dir_base); 