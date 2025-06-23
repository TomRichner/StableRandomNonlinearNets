% consolidate_parameter_space_results.m
%
% This script is designed to be run after a 'network_parameter_space_analysis'
% simulation has been stopped prematurely, leaving behind a 'temp_results'
% directory with unconsolidated batch files. It will load these batch files,
% assemble them into a final results structure, and save it in the format
% expected by downstream analysis scripts.

close all;
clear all;
clc;

%% Get the base directory for the analysis
fprintf('Please select the main experiment directory to consolidate...\n');
base_dir = uigetdir(pwd, 'Select Experiment Directory');
if isequal(base_dir, 0)
    fprintf('User cancelled. Exiting.\n');
    return;
end
fprintf('Selected directory: %s\n', base_dir);

%% Find the temp_results directory
temp_output_dir = fullfile(base_dir, 'temp_results');
if ~exist(temp_output_dir, 'dir')
    fprintf('Could not find "temp_results" directory in the selected folder.\n');
    fprintf('Please select the top-level experiment directory (e.g., "space_search_...").\n');
    return;
end

%% Find all batch files and infer experiment parameters
batch_files = dir(fullfile(temp_output_dir, 'batch_results_*.mat'));
if isempty(batch_files)
    fprintf('No batch files found in "temp_results". Nothing to consolidate.\n');
    return;
end
fprintf('Found %d batch files to process.\n\n', length(batch_files));

% --- Load the first batch file to learn about the conditions ---
try
    first_batch_data = load(fullfile(batch_files(1).folder, batch_files(1).name));
    condition_names = fieldnames(first_batch_data.batch_results_by_condition);
    num_conditions = length(condition_names);
    fprintf('Inferred %d conditions from batch file: %s\n', num_conditions, strjoin(condition_names, ', '));
catch
    fprintf('Error: Could not read batch file format. The file may be corrupt.\n');
    return;
end

% --- Determine total number of combinations from all batch files ---
batch_size = 500; % Must match the batch_size in the simulation script
max_batch_num = 0;
num_in_last_batch = 0;

for i = 1:length(batch_files)
    [~, name, ~] = fileparts(batch_files(i).name);
    num = sscanf(name, 'batch_results_%d');
    if num > max_batch_num
        max_batch_num = num;
        loaded_data = load(fullfile(temp_output_dir, batch_files(i).name));
        % Get size from the first condition's results, assuming all are the same
        num_in_last_batch = length(loaded_data.batch_results_by_condition.(condition_names{1}));
    end
end

if max_batch_num == 0
    fprintf('Could not determine batch numbering. Cannot proceed.\n');
    return;
end

num_combinations = (max_batch_num - 1) * batch_size + num_in_last_batch;
fprintf('Inferred %d total combinations from %d batch files.\n', num_combinations, length(batch_files));

%% Loop through each condition and consolidate its data
for c_idx = 1:num_conditions
    condition_name = condition_names{c_idx};
    output_dir = fullfile(base_dir, condition_name);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    fprintf('======================================================\n');
    fprintf('===== Consolidating Condition: %s =====\n', upper(condition_name));
    
    final_save_filename = fullfile(output_dir, sprintf('param_space_results_%s.mat', condition_name));
    if exist(final_save_filename, 'file')
        choice = questdlg(sprintf('Final results file already exists for "%s". Overwrite?', condition_name), ...
            'File Exists', 'Yes', 'No', 'No');
        if strcmp(choice, 'No')
            fprintf('Skipping condition "%s" as per user request.\n\n', condition_name);
            continue;
        end
    end

    % --- Consolidate results for this condition using a matfile object---
    fprintf('Consolidating data for "%s"...\n', condition_name);
    matObj = matfile(final_save_filename, 'Writable', true);
    matObj.results = cell(num_combinations, 1); % pre-allocate
    all_batches_found_for_cond = true;

    for i = 1:length(batch_files)
        batch_filename = fullfile(batch_files(i).folder, batch_files(i).name);
        [~, name, ~] = fileparts(batch_filename);
        batch_idx = sscanf(name, 'batch_results_%d');

        try
            loaded_data = load(batch_filename);
            batch_slice = loaded_data.batch_results_by_condition.(condition_name);
            
            start_idx = (batch_idx - 1) * batch_size + 1;
            end_idx = start_idx + length(batch_slice) - 1;
            
            matObj.results(start_idx:end_idx, 1) = batch_slice;
        catch
            fprintf('Warning: Could not read or place data for "%s" from batch file: %s\n', condition_name, batch_files(i).name);
            all_batches_found_for_cond = false;
        end
    end
    
    % --- Save final file with partial metadata ---
    fprintf('Saving final file with metadata...\n');
    metadata = struct();
    metadata.NOTE = 'This file was created by consolidate_parameter_space_results.m from an interrupted run.';
    metadata.consolidation_date = datestr(now);
    metadata.num_combinations = num_combinations;
    metadata.original_run_directory = base_dir;
    
    matObj.metadata = metadata;
    fprintf('Successfully saved consolidated results to:\n%s\n\n', final_save_filename);
end

%% Final cleanup
fprintf('======================================================\n');
fprintf('===== Consolidation Complete for all conditions =====\n');

choice = questdlg('Delete the temporary results folder?', 'Clean Up', 'Yes', 'No', 'No');
if strcmp(choice, 'Yes')
    rmdir(temp_output_dir, 's');
    fprintf('Temporary results directory has been deleted.\n');
else
    fprintf('Temporary results directory was not deleted.\n');
end

fprintf('======================================================\n');
beep; 