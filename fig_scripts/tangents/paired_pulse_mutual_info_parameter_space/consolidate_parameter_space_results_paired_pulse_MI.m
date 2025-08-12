% consolidate_parameter_space_results_paired_pulse_MI.m
%
% Consolidate 'temp_results' batch files from an interrupted paired-pulse
% parameter space exploration into final per-condition .mat files that match
% the format expected by plotting scripts.

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
    fprintf('Please select the top-level experiment directory (e.g., "paired_pulse_...").\n');
    return;
end

%% Find all batch files and infer experiment parameters
batch_files = dir(fullfile(temp_output_dir, 'batch_results_*.mat'));
if isempty(batch_files)
    fprintf('No batch files found in "temp_results". Nothing to consolidate.\n');
    return;
end
fprintf('Found %d batch files to process.\n\n', length(batch_files));

% Load the first batch to infer conditions
try
    first_batch_data = load(fullfile(batch_files(1).folder, batch_files(1).name));
    condition_names = fieldnames(first_batch_data.batch_results_by_condition);
    num_conditions = length(condition_names);
    fprintf('Inferred %d conditions from batch file: %s\n', num_conditions, strjoin(condition_names, ', '));
catch
    fprintf('Error: Could not read batch file format. The file may be corrupt.\n');
    return;
end

% Determine total number of combinations
prompt = {'Enter the batch size used in the original simulation:'};
dlgtitle = 'Input Batch Size';
dims = [1 50];
definput = {'25'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    fprintf('User cancelled. Exiting.\n');
    return;
end
batch_size = str2double(answer{1});
if isnan(batch_size) || batch_size <= 0 || floor(batch_size) ~= batch_size
    fprintf('Invalid batch size entered. Must be a positive integer. Exiting.\n');
    return;
end

batch_numbers = zeros(length(batch_files), 1);
for i = 1:length(batch_files)
    [~, name, ~] = fileparts(batch_files(i).name);
    scanned_num = sscanf(name, 'batch_results_%d');
    if ~isempty(scanned_num)
        batch_numbers(i) = scanned_num(1);
    end
end

if all(batch_numbers == 0)
    fprintf('Could not determine batch numbering from filenames. Cannot proceed.\n');
    return;
end

[max_batch_num, max_idx] = max(batch_numbers); %#ok<ASGLU>
last_batch_data = load(fullfile(batch_files(max_idx).folder, batch_files(max_idx).name));
num_in_last_batch = length(last_batch_data.batch_results_by_condition.(condition_names{1}));

num_combinations = (max(batch_numbers) - 1) * batch_size + num_in_last_batch;
fprintf('Inferred %d total combinations from %d batch files (using batch size of %d).\n\n', num_combinations, length(batch_files), batch_size);

%% Pre-allocate results in memory for all conditions
fprintf('======================================================\n');
fprintf('===== Allocating Memory For In-Memory Consolidation =====\n');
fprintf('WARNING: This will load all results into RAM. This may fail if the dataset is very large.\n');

results_in_memory = struct();
for c_idx = 1:num_conditions
    condition_name = condition_names{c_idx};
    results_in_memory.(condition_name) = cell(num_combinations, 1);
    fprintf('Allocated memory for "%s" (%d entries).\n', condition_name, num_combinations);
end
fprintf('======================================================\n\n');

%% Load all batch file data into memory
fprintf('===== Loading all batch files into RAM =====\n');
all_batches_found_and_loaded = true;
for i = 1:length(batch_files)
    batch_filename = fullfile(batch_files(i).folder, batch_files(i).name);
    batch_idx = batch_numbers(i);

    if batch_idx == 0
        fprintf('Skipping file with unrecognized batch number: %s\n', batch_files(i).name);
        continue;
    end

    fprintf('Reading batch file %d of %d (%s)...\n', i, length(batch_files), batch_files(i).name);

    try
        loaded_data = load(batch_filename, 'batch_results_by_condition');
        start_idx = (batch_idx - 1) * batch_size + 1;

        for c_idx = 1:num_conditions
            condition_name = condition_names{c_idx};
            if isfield(loaded_data.batch_results_by_condition, condition_name)
                batch_slice = loaded_data.batch_results_by_condition.(condition_name);
                end_idx = start_idx + length(batch_slice) - 1;

                if end_idx <= num_combinations
                    results_in_memory.(condition_name)(start_idx:end_idx, 1) = batch_slice;
                else
                    fprintf('  Warning: Index out of bounds for condition "%s". Skipping placement.\n', condition_name);
                end
            end
        end
    catch ME
        fprintf('  Error reading or processing file: %s. Error: %s\n', batch_files(i).name, ME.message);
        all_batches_found_and_loaded = false;
    end
end
fprintf('======================================================\n\n');

%% Write final consolidated files from memory
fprintf('===== Writing consolidated files to disk =====\n');
for c_idx = 1:num_conditions
    condition_name = condition_names{c_idx};
    output_dir = fullfile(base_dir, condition_name);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    final_save_filename = fullfile(output_dir, sprintf('param_space_results_%s.mat', condition_name));

    choice = 'Yes';
    if exist(final_save_filename, 'file')
        choice = questdlg(sprintf('Final results file already exists for "%s". Overwrite?', condition_name), ...
            'File Exists', 'Yes', 'No', 'No');
    end

    if strcmp(choice, 'Yes')
        fprintf('Writing final file for "%s"...\n', condition_name);

        results = results_in_memory.(condition_name); %#ok<NASGU>

        metadata = struct();
        metadata.NOTE = 'This file was created by consolidate_parameter_space_results_paired_pulse_MI.m from an interrupted run.';
        metadata.consolidation_date = datestr(now);
        metadata.num_combinations = num_combinations;
        metadata.original_run_directory = base_dir;

        vars_to_save = {'results', 'metadata'};

        try
            save(final_save_filename, vars_to_save{:}, '-v7.3');
            fprintf('Successfully saved consolidated results to:\n%s\n\n', final_save_filename);
        catch ME
            fprintf('ERROR: Could not save the final file for "%s".\n', condition_name);
            fprintf('The error was: %s\n', ME.message);
            fprintf('The variable might be too large to save. You may need to use a matfile-based approach for consolidation.\n\n');
        end
    else
        fprintf('Skipping condition "%s" as per user request.\n\n', condition_name);
    end
end

%% Final cleanup
fprintf('======================================================\n');
fprintf('===== Consolidation Complete =====\n');

choice = questdlg('Delete the temporary results folder?', 'Clean Up', 'Yes', 'No', 'No');
if strcmp(choice, 'Yes') && all_batches_found_and_loaded
    rmdir(temp_output_dir, 's');
    fprintf('Temporary results directory has been deleted.\n');
elseif strcmp(choice, 'Yes') && ~all_batches_found_and_loaded
    fprintf('Did not delete temp folder because some batches were missing or failed to load.\n');
else
    fprintf('Temporary results directory was not deleted.\n');
end

fprintf('======================================================\n');
beep;


