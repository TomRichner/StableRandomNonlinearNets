% Script for plotting sensitivity analysis results
% This script loads the sensitivity analysis results and creates visualizations
% of the largest Lyapunov exponent (LLE) across parameter ranges

clear;
clc;
close all;

% Add the root directory to the path to access model functions
addpath(fullfile(fileparts(mfilename('fullpath')),'../..'));

%% Set up paths and parameters
sensitivity_results_dir = fullfile(fileparts(mfilename('fullpath')), '../../sensitivity_results');
output_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Define LLE histogram parameters
lle_range = [-0.25, 0.25];
n_bins = 30;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins), inf];

%% Load summary file to get list of analyzed parameters
summary_file = fullfile(sensitivity_results_dir, 'sensitivity_analysis_summary.mat');
if ~exist(summary_file, 'file')
    error('Summary file not found: %s', summary_file);
end

fprintf('Loading summary file: %s\n', summary_file);
summary_data = load(summary_file);

% Get parameter names from the summary
if isfield(summary_data, 'summary_data') && isfield(summary_data.summary_data, 'all_results_summary')
    param_names = fieldnames(summary_data.summary_data.all_results_summary);
else
    % Fallback: search for sensitivity_*.mat files in the directory
    files = dir(fullfile(sensitivity_results_dir, 'sensitivity_*.mat'));
    param_names = {};
    for i = 1:length(files)
        if ~strcmp(files(i).name, 'sensitivity_analysis_summary.mat')
            param_name = strrep(files(i).name, 'sensitivity_', '');
            param_name = strrep(param_name, '.mat', '');
            param_names{end+1} = param_name;
        end
    end
end

fprintf('Found %d parameters to plot: %s\n', length(param_names), strjoin(param_names, ', '));

%% Process each parameter
for i = 1:length(param_names)
    param_name = param_names{i};
    fprintf('\n--- Processing %s (%d/%d) ---\n', param_name, i, length(param_names));
    
    % Load parameter-specific results
    param_file = fullfile(sensitivity_results_dir, sprintf('sensitivity_%s.mat', param_name));
    if ~exist(param_file, 'file')
        fprintf('Warning: File not found: %s\n', param_file);
        continue;
    end
    
    try
        % Call the plotting function
        sensitivity_plot(param_file, lle_bins, output_dir);
        fprintf('Successfully plotted %s\n', param_name);
    catch ME
        fprintf('Error plotting %s: %s\n', param_name, ME.message);
        if ~isempty(ME.stack)
            fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

fprintf('\n=== Sensitivity plotting complete ===\n');
fprintf('Plots saved to: %s\n', output_dir); 