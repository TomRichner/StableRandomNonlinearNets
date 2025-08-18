% network_parameter_space_analysis_paired_pulse_MI.m
% Explore SRNN parameter space with paired-pulse stimulus and compute MI vs delay.
% Mirrors fig_scripts/SRNN_paper/parameter_space_analysis/network_parameter_space_analysis.m
% but uses SRNN_caller_wrapped_for_paired_pulse_MI for stimulus and MI.

close all;
clear all;
clc;

%% Analysis Conditions to Compare
conditions = { ...
    struct('name', 'no_adaptation', 'n_a_E_val', 0, 'n_b_E_val', 0), ...
    struct('name', 'sfa_only',      'n_a_E_val', 1, 'n_b_E_val', 0), ...
    struct('name', 'std_only',      'n_a_E_val', 0, 'n_b_E_val', 1), ...
    struct('name', 'sfa_and_std',   'n_a_E_val', 1, 'n_b_E_val', 1) ...
};

%% Parameters for Grid Search
params_for_grid = {'EE_factor', 'mean_weight', 'tau_a_E_2', 'tau_b_E_2'};

%% Analysis Parameters
n_levels = 5; % Number of values to test for each parameter in the grid
note = 'ppMI_param_space';
% Timestamp for folder name
dt_str = lower(strrep(datestr(now, 'mmm_dd_yy_hh_MM_AM'), ':', '_'));

%% Create an ABSOLUTE base directory and be sure it exists
output_dir_base = fullfile(pwd, ['paired_pulse_' note '_nLevs_' num2str(n_levels) '_' dt_str]);
if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

%% copy current mfile (this file) into it
copyfile([mfilename('fullpath') '.m'], output_dir_base);

%% Loop through conditions
overall_super_start_time = tic;
all_conditions_summary = struct();


fprintf('\n\n======================================================\n');
fprintf('===== Running Paired-Pulse MI Analysis for All Conditions =====\n');
fprintf('======================================================\n\n');

%% Default Simulation Parameters
p_default.fs = 500;
p_default.n = 10;
p_default.EE_factor = 1.0;
p_default.IE_factor = 1.0;
p_default.EI = 0.7;
p_default.E_self = 0.0;
p_default.mean_weight = 0.5;
p_default.DC = 0.1;
p_default.mean_in_out_degree = 5;
p_default.tau_a_E_2 = 15;
p_default.tau_b_E_2 = 2;
p_default.tau_STD = 0.5;
p_default.c_SFA_factor = 0.5;

%% Parameter Ranges for Grid Search
ranges.fs = [250, 2000];
ranges.n = [10, 200];
ranges.EE_factor = [0.5 1.5];
ranges.IE_factor = [0.05, 2];
ranges.EI = [1/p_default.n, 1.0];
ranges.E_self = [0.0, 0.5];
ranges.mean_weight = [0.25 0.75];
ranges.DC = [0, 4];
ranges.mean_in_out_degree = [1.5, p_default.n-1];
ranges.tau_a_E_2 = [1, 15];
ranges.tau_b_E_2 = [1, 15];
ranges.tau_STD = [0.2, 1];
ranges.c_SFA_factor = [0.0, 4.0];

% Validate that all specified parameters for the grid exist in ranges
all_possible_params = fieldnames(ranges);
invalid_params = setdiff(params_for_grid, all_possible_params);
if ~isempty(invalid_params)
    error('Invalid parameter names in params_for_grid: %s', strjoin(invalid_params, ', '));
end

% Create a temporary directory for intermediate results at the top level
temp_output_dir = fullfile(output_dir_base, 'temp_results');
if ~exist(temp_output_dir, 'dir')
    mkdir(temp_output_dir);
end

%% Generate Parameter Grid
param_vectors = cell(1, length(params_for_grid));
for i = 1:length(params_for_grid)
    param_name = params_for_grid{i};
    param_range = ranges.(param_name);
    if ismember(param_name, {'n', 'mean_in_out_degree'})
        param_vectors{i} = round(linspace(param_range(1), param_range(2), n_levels));
    else
        param_vectors{i} = linspace(param_range(1), param_range(2), n_levels);
    end
end

grid_cells = cell(size(param_vectors));
[grid_cells{:}] = ndgrid(param_vectors{:});

num_combinations = numel(grid_cells{1});
all_configs = cell(num_combinations, 1);
for i = 1:num_combinations
    config = struct();
    for j = 1:length(params_for_grid)
        config.(params_for_grid{j}) = grid_cells{j}(i);
    end
    all_configs{i} = config;
end

% Randomize the order of execution by shuffling the configs
shuffled_configs = all_configs(randperm(num_combinations));

fprintf('=== Paired-Pulse Parameter Space Exploration Started for All Conditions ===\n');
fprintf('Grid parameters: %s\n', strjoin(params_for_grid, ', '));
fprintf('Levels per parameter: %d\n', n_levels);
fprintf('Total combinations to simulate: %d\n', num_combinations);
fprintf('=========================================================\n\n');

%% Run Simulations in Batches, Interleaving Conditions
batch_size = 10; % Number of parameter configs per batch
num_batches = ceil(num_combinations / batch_size);
num_conditions = length(conditions);
fprintf('Running %d combinations for %d conditions in %d batches.\n', num_combinations, num_conditions, num_batches);

for batch_idx = 1:num_batches
    batch_save_filename = fullfile(temp_output_dir, sprintf('batch_results_%d.mat', batch_idx));
    if exist(batch_save_filename, 'file')
        fprintf('Batch %d/%d already completed. Skipping.\n', batch_idx, num_batches);
        continue;
    end

    start_idx = (batch_idx - 1) * batch_size + 1;
    end_idx = min(batch_idx * batch_size, num_combinations);

    configs_in_batch = shuffled_configs(start_idx:end_idx);
    current_batch_size = length(configs_in_batch);
    num_jobs_in_batch = current_batch_size * num_conditions;

    fprintf('--- Running Batch %d/%d (Configs %d-%d for all %d conditions) ---\n', batch_idx, num_batches, start_idx, end_idx, num_conditions);

    % --- Create a list of all jobs for the parfor loop ---
    jobs_to_run = cell(num_jobs_in_batch, 1);
    job_counter = 1;
    for config_k = 1:current_batch_size
        for cond_k = 1:num_conditions
            job = struct();
            job.config = configs_in_batch{config_k};
            job.condition = conditions{cond_k};
            job.config_run_idx = start_idx + config_k - 1;
            job.condition_idx = cond_k;
            jobs_to_run{job_counter} = job;
            job_counter = job_counter + 1;
        end
    end

    % --- Run all jobs for the batch in parallel ---
    parallel_results = cell(num_jobs_in_batch, 1);
    parfor j = 1:num_jobs_in_batch
    % for j = 1:num_jobs_in_batch
        current_job = jobs_to_run{j};
        current_p = p_default;

        % Set condition-specific parameters
        current_p.n_a_E = current_job.condition.n_a_E_val;
        current_p.n_b_E = current_job.condition.n_b_E_val;

        % Set grid parameters
        config_fields = fieldnames(current_job.config);
        for f_idx = 1:length(config_fields)
            field = config_fields{f_idx};
            current_p.(field) = current_job.config.(field);
        end

        sim_seed = current_job.config_run_idx + (current_job.condition_idx-1)*num_combinations; % Globally unique seed

        % try
            result_sim = SRNN_caller_wrapped_for_paired_pulse_MI(...
                sim_seed, current_p.n, current_p.EE_factor, current_p.IE_factor, ...
                current_p.EI, current_p.E_self, current_p.mean_weight, current_p.DC, ...
                current_p.mean_in_out_degree, current_p.tau_a_E_2, current_p.tau_b_E_2, ...
                current_p.tau_STD, current_p.c_SFA_factor, current_p.n_a_E, ...
                current_p.n_b_E, current_p.fs);

            result_sim.success = true;

            result_to_save = struct();
            result_to_save.simulation_output = result_sim;
            result_to_save.parameters = current_p;
            result_to_save.config = current_job.config;
            result_to_save.success = true;
            result_to_save.seed = sim_seed;

            parallel_results{j} = result_to_save;

        % catch ME
        %     result_to_save = struct(...
        %         'success', false, 'error', ME, 'error_message', ME.message, ...
        %         'parameters', current_p, 'config', current_job.config, 'seed', sim_seed);
        % 
        %     parallel_results{j} = result_to_save;
        % end
    end

    % --- Organize flat results back into per-condition structs ---
    batch_results_by_condition = struct();
    for cond_k = 1:num_conditions
        batch_results_by_condition.(conditions{cond_k}.name) = cell(current_batch_size, 1);
    end

    for j = 1:num_jobs_in_batch
        condition_name = jobs_to_run{j}.condition.name;
        config_idx_in_batch = floor((j-1)/num_conditions) + 1;
        batch_results_by_condition.(condition_name){config_idx_in_batch} = parallel_results{j};
    end

    fprintf('--- Batch %d/%d finished. Saving results... ---\n', batch_idx, num_batches);
    save(batch_save_filename, 'batch_results_by_condition', '-v7.3');
    fprintf('Batch results saved to %s\n\n', batch_save_filename);
end

%% Consolidate results and Finalize
fprintf('--- Consolidating all batch results ---\n');

% Create final output directories and matfile objects for each condition
matObj = struct();
for c2_idx = 1:num_conditions
    condition_name2 = conditions{c2_idx}.name;
    output_dir2 = fullfile(output_dir_base, condition_name2);
    if ~exist(output_dir2, 'dir'), mkdir(output_dir2); end

    save_filename2 = fullfile(output_dir2, sprintf('param_space_results_%s.mat', condition_name2));
    if exist(save_filename2, 'file'), delete(save_filename2); end

    matObj.(condition_name2) = matfile(save_filename2, 'Writable', true);
    matObj.(condition_name2).results = cell(num_combinations, 1);
end

all_batches_found = true;
for batch_idx = 1:num_batches
    batch_save_filename = fullfile(temp_output_dir, sprintf('batch_results_%d.mat', batch_idx));
    if exist(batch_save_filename, 'file')
        loaded_data = load(batch_save_filename);

        start_idx = (batch_idx - 1) * batch_size + 1;
        end_idx = min(batch_idx * batch_size, num_combinations);

        if isfield(loaded_data, 'batch_results_by_condition')
            for c2_idx = 1:num_conditions
                condition_name2 = conditions{c2_idx}.name;
                batch_slice = loaded_data.batch_results_by_condition.(condition_name2);
                matObj.(condition_name2).results(start_idx:end_idx, 1) = batch_slice;
            end
        else
            fprintf('Warning: Corrupt batch file for batch %d. Skipping.\n', batch_idx);
            all_batches_found = false;
        end
    else
        fprintf('Warning: Result file for batch %d not found!\n', batch_idx);
        all_batches_found = false;
    end
end

if all_batches_found
    fprintf('All batch result files consolidated successfully.\n');
else
    fprintf('Warning: Some batch result files were missing or corrupt.\n');
end

% Save metadata into each final MAT-file and store summary
total_elapsed_time = toc(overall_super_start_time); % Calculate total time after batches
all_conditions_summary = struct();
for c2_idx = 1:num_conditions
    condition_name2 = conditions{c2_idx}.name;

    metadata = struct();
    metadata.params_for_grid = params_for_grid;
    metadata.param_vectors = param_vectors;
    metadata.shuffled_configs = shuffled_configs;
    
    % Create a clean p_default for each condition's metadata
    p_default_for_meta = p_default;
    p_default_for_meta.n_a_E = conditions{c2_idx}.n_a_E_val;
    p_default_for_meta.n_b_E = conditions{c2_idx}.n_b_E_val;
    metadata.p_default = p_default_for_meta;

    metadata.n_levels = n_levels;
    metadata.num_combinations = num_combinations;
    metadata.elapsed_time_sec = total_elapsed_time;
    metadata.analysis_date = datestr(now);
    metadata.condition = conditions{c2_idx};

    matObj.(condition_name2).metadata = metadata;

    summary = struct();
    summary.metadata = metadata;
    summary.save_filename = matObj.(condition_name2).Properties.Source;
    all_conditions_summary.(condition_name2) = summary;
end

fprintf('Final results and metadata saved for all conditions.\n');

% Clean up temporary directory only if consolidation was successful
if all_batches_found
    rmdir(temp_output_dir, 's');
    fprintf('Temporary results directory cleaned up.\n');
else
    fprintf('Temporary results directory retained due to missing batches.\n');
end

%% Final summary
total_elapsed = toc(overall_super_start_time);
fprintf('=== OVERALL PAIRED-PULSE PARAMETER SPACE EXPLORATION COMPLETE ===\n');
fprintf('Total duration: %.2f hours\n', total_elapsed/3600);
fprintf('Conditions analyzed: %d\n', length(conditions));
fprintf('==================================================\n');

% Save final summary
final_summary_filename = fullfile(output_dir_base, 'parameter_space_analysis_summary.mat');
save(final_summary_filename, 'all_conditions_summary', 'total_elapsed', 'conditions');
fprintf('Final summary of all conditions saved to %s\n', final_summary_filename);
fprintf('==================================================\n');

beep;pause(0.5);beep;pause(0.2);beep % done, wake up


