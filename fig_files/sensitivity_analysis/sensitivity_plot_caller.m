% Script for plotting sensitivity analysis results
% This script loads the sensitivity analysis results and creates visualizations
% of the largest Lyapunov exponent (LLE) across parameter ranges

clear;
clc;
close all;


%% Set up paths and parameters
sensitivity_results_dir = '/Users/tom/Desktop/local_code/adaptation_paper/dual_adaptation_random_matrix_theory/fig_files/sensitivity_analysis/sensitivity_results';
output_dir_base = fullfile(fileparts(mfilename('fullpath')), 'plots');

% Create base output directory if it doesn't exist
if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

% Define LLE histogram parameters
lle_range = [-0.25, 0.25];
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

%% Process each condition
for c_idx = 1:length(conditions)
    current_condition = conditions{c_idx};
    condition_name = current_condition.name;
    
    fprintf('\n\n--- Processing Condition: %s ---\n', upper(condition_name));

    condition_results_dir = fullfile(sensitivity_results_dir, condition_name);
    output_dir = fullfile(output_dir_base, condition_name);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    % Get parameter names from the summary for this condition
    param_names = {};
    if isfield(summary_data, 'all_conditions_summary') && isfield(summary_data.all_conditions_summary, condition_name)
        param_names = fieldnames(summary_data.all_conditions_summary.(condition_name));
    else
        fprintf('Warning: Could not find summary for condition "%s" in summary file. Skipping.\n', condition_name);
        continue;
    end
    
    if isempty(param_names)
        fprintf('No parameters found for condition "%s". Skipping.\n', condition_name);
        continue;
    end

    fprintf('Found %d parameters to plot for %s: %s\n', length(param_names), condition_name, strjoin(param_names, ', '));

    %% Process each parameter for the current condition
    for i = 1:length(param_names)
        param_name = param_names{i};
        fprintf('\n--- Processing %s (%d/%d) ---\n', param_name, i, length(param_names));
        
        % Load parameter-specific results
        param_file = fullfile(condition_results_dir, sprintf('sensitivity_%s.mat', param_name));
        if ~exist(param_file, 'file')
            fprintf('Warning: File not found: %s\n', param_file);
            continue;
        end
        
        try
            % Call the plotting function
            sensitivity_plot(param_file, lle_bins, output_dir);
            fprintf('Successfully plotted %s for condition %s\n', param_name, condition_name);
        catch ME
            fprintf('Error plotting %s for condition %s: %s\n', param_name, condition_name, ME.message);
            if ~isempty(ME.stack)
                fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
    end
end

fprintf('\n=== Sensitivity plotting complete ===\n');
fprintf('Plots saved to: %s\n', output_dir_base); 