% Script for exploring the parameter space of the SRNN model.
% This script generates a grid of parameter combinations, simulates each point
% in a random order, and saves the results. This is intended for running
% long explorations of the model's behavior across different parameter regimes.

close all;
clear all;
clc;

%% Analysis Conditions to Compare
conditions = { ...
    struct('name', 'no_adaptation', 'n_a_E_val', 0, 'n_b_E_val', 0), ...
    struct('name', 'sfa_only',      'n_a_E_val', 3, 'n_b_E_val', 0), ...
    struct('name', 'std_only',      'n_a_E_val', 0, 'n_b_E_val', 2), ...
    struct('name', 'sfa_and_std',   'n_a_E_val', 3, 'n_b_E_val', 2) ...
};

%% Parameters for Grid Search
% Specify which parameters to vary in the grid search.
% The script will generate all combinations of the specified levels for these parameters.
% Parameters not in this list will be held at their default values.
params_for_grid = {'EE_factor', 'mean_weight', 'EI', 'IE_factor','mean_in_out_degree','tau_a_E_2'};

%% Analysis Parameters
n_levels = 2; % Number of values to test for each parameter in the grid
note = 'param_space_exploration_LLE_SR_EE_W_n';

% Timestamp for folder name
dt_str = lower(strrep(datestr(now, 'mmm_dd_yy_hh_MM_AM'), ':', '_'));

%% Create an ABSOLUTE base directory and be sure it exists
output_dir_base = fullfile(pwd, ...
    ['space_search_' note '_nLevs_' num2str(n_levels) '_' dt_str]);

if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

%% copy current mfile (this file) into it
copyfile([mfilename('fullpath') '.m'], output_dir_base);

%%  loop through conditions
overall_super_start_time = tic;
all_conditions_summary = struct();

for c_idx = 1:length(conditions)
    current_condition = conditions{c_idx};
    fprintf('\n\n======================================================\n');
    fprintf('===== Running Analysis for Condition: %s =====\n', upper(current_condition.name));
    fprintf('===== n_a_E = %d, n_b_E = %d =====\n', current_condition.n_a_E_val, current_condition.n_b_E_val);
    fprintf('======================================================\n\n');
    
    %% Default Simulation Parameters
    p_default.fs = 1000;
    p_default.n = 10;
    p_default.EE_factor = 1.0;
    p_default.IE_factor = 1.0;
    p_default.EI = 0.7;
    p_default.E_self = 0.0;
    p_default.mean_weight = 0.5;
    p_default.DC = 0.1;
    p_default.mean_in_out_degree = 5;
    p_default.tau_a_E_2 = 15;
    p_default.tau_b_E_2 = 9;
    p_default.tau_STD = 0.5;
    p_default.c_SFA_factor = 0.5;
    p_default.n_a_E = current_condition.n_a_E_val;
    p_default.n_b_E = current_condition.n_b_E_val;
    
    %% Parameter Ranges for Grid Search
    ranges.fs = [250, 2000];
    ranges.n = [10, 200];
    ranges.EE_factor = [0.05, 4];
    ranges.IE_factor = [0.05, 4];
    ranges.EI = [1/p_default.n, 1.0];
    ranges.E_self = [0.0, 0.5];
    ranges.mean_weight = [1/p_default.n, 2];
    ranges.DC = [0, 4];
    ranges.mean_in_out_degree = [1.5, p_default.n-1];
    ranges.tau_a_E_2 = [6, 30];
    ranges.tau_b_E_2 = [2, 30];
    ranges.tau_STD = [0.2, 1];
    ranges.c_SFA_factor = [0.0, 4.0];
    
    % Validate that all specified parameters for the grid exist in ranges
    all_possible_params = fieldnames(ranges);
    invalid_params = setdiff(params_for_grid, all_possible_params);
    if ~isempty(invalid_params)
        error('Invalid parameter names in params_for_grid: %s', strjoin(invalid_params, ', '));
    end
    
    output_dir = fullfile(output_dir_base, current_condition.name);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Create a temporary directory for intermediate results
    temp_output_dir = fullfile(output_dir, 'temp_results');
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
    
    fprintf('=== SRNN Parameter Space Exploration Started for %s ===\n', current_condition.name);
    fprintf('Grid parameters: %s\n', strjoin(params_for_grid, ', '));
    fprintf('Levels per parameter: %d\n', n_levels);
    fprintf('Total combinations to simulate: %d\n', num_combinations);
    fprintf('=========================================================\n\n');
    
    %% Run Simulations in Batches
    condition_start_time = tic;
    batch_size = 500;
    num_batches = ceil(num_combinations / batch_size);
    fprintf('Running %d combinations in %d batches of up to size %d.\n', num_combinations, num_batches, batch_size);

    for batch_idx = 1:num_batches
        batch_save_filename = fullfile(temp_output_dir, sprintf('batch_results_%d.mat', batch_idx));
        if exist(batch_save_filename, 'file')
            fprintf('Batch %d/%d already completed. Skipping.\n', batch_idx, num_batches);
            continue;
        end

        start_idx = (batch_idx - 1) * batch_size + 1;
        end_idx = min(batch_idx * batch_size, num_combinations);
        current_batch_indices = start_idx:end_idx;
        current_batch_size = length(current_batch_indices);

        fprintf('--- Running Batch %d/%d (Simulations %d to %d) ---\n', batch_idx, num_batches, start_idx, end_idx);
        
        batch_results = cell(current_batch_size, 1);
        
        parfor k = 1:current_batch_size
            config_idx = current_batch_indices(k);
            current_config = shuffled_configs{config_idx};

            current_p = p_default;
            config_fields = fieldnames(current_config);
            for f_idx = 1:length(config_fields)
                field = config_fields{f_idx};
                current_p.(field) = current_config.(field);
            end
            
            sim_seed = config_idx; % Use a unique seed for each config across all batches
            
            run_start = tic;
            try
                result_sim = SRNN_caller_wrapped_for_sensitivity_dual_stage_rndWalk(...
                    sim_seed, ...
                    current_p.n, ...
                    current_p.EE_factor, ...
                    current_p.IE_factor, ...
                    current_p.EI, ...
                    current_p.E_self, ...
                    current_p.mean_weight, ...
                    current_p.DC, ...
                    current_p.mean_in_out_degree, ...
                    current_p.tau_a_E_2, ...
                    current_p.tau_b_E_2, ...
                    current_p.tau_STD, ...
                    current_p.c_SFA_factor, ...
                    current_p.n_a_E, ...
                    current_p.n_b_E, ...
                    current_p.fs);
                
                run_duration = toc(run_start);
                result_sim.run_duration = run_duration;
                result_sim.success = true;
                
                result_to_save = struct();
                result_to_save.simulation_output = result_sim;
                result_to_save.parameters = current_p;
                result_to_save.config = current_config;
                result_to_save.success = true;
                result_to_save.seed = sim_seed;
                result_to_save.run_order_index = config_idx;
                
                batch_results{k} = result_to_save;
                
            catch ME
                run_duration = toc(run_start);
                result_to_save = struct(...
                    'success', false, ...
                    'error', ME, ...
                    'error_message', ME.message, ...
                    'parameters', current_p, ...
                    'config', current_config, ...
                    'seed', sim_seed, ...
                    'run_duration', run_duration, ...
                    'run_order_index', config_idx);
                
                batch_results{k} = result_to_save;
            end
        end
        
        fprintf('--- Batch %d/%d finished. Saving results... ---\n', batch_idx, num_batches);
        save(batch_save_filename, 'batch_results', '-v7.3');
        fprintf('Batch results saved to %s\n', batch_save_filename);
    end
    
    condition_elapsed_time = toc(condition_start_time);
    
    %% Consolidate results and Finalize
    fprintf('--- Consolidating results for condition: %s ---\n', current_condition.name);

    % Create / open a matfile for streaming writes
    save_filename = fullfile(output_dir, sprintf('param_space_results_%s.mat', current_condition.name));
    if exist(save_filename, 'file')
        delete(save_filename); % start fresh to avoid stale content
    end
    matObj = matfile(save_filename, 'Writable', true);
    matObj.results = cell(num_combinations, 1); % pre-allocate on disk

    all_batches_found = true;
    for batch_idx = 1:num_batches
        batch_save_filename = fullfile(temp_output_dir, sprintf('batch_results_%d.mat', batch_idx));
        if exist(batch_save_filename, 'file')
            loaded_data = load(batch_save_filename);

            start_idx = (batch_idx - 1) * batch_size + 1;
            end_idx = min(batch_idx * batch_size, num_combinations);

            if isfield(loaded_data, 'batch_results')
                % Stream-write this batch slice directly to the MAT-file
                matObj.results(start_idx:end_idx, 1) = loaded_data.batch_results;
            else
                fprintf('Warning: Corrupt or old batch file format for batch %d. Skipping.\n', batch_idx);
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

    % Save metadata into the same MAT-file
    metadata = struct();
    metadata.params_for_grid = params_for_grid;
    metadata.param_vectors = param_vectors;
    metadata.shuffled_configs = shuffled_configs;
    metadata.p_default = p_default;
    metadata.n_levels = n_levels;
    metadata.num_combinations = num_combinations;
    metadata.elapsed_time_sec = condition_elapsed_time;
    metadata.analysis_date = datestr(now);
    metadata.condition = current_condition;

    matObj.metadata = metadata;

    fprintf('\n--- Analysis for condition %s finished in %.2f hours ---\n', ...
        current_condition.name, condition_elapsed_time/3600);
    fprintf('Results saved to %s\n\n', save_filename);

    % Clean up temporary directory only if consolidation was successful
    if all_batches_found
        rmdir(temp_output_dir, 's');
        fprintf('Temporary results directory cleaned up.\n');
    else
        fprintf('Temporary results directory retained due to missing batches.\n');
    end

    % Store summary for this condition
    summary = struct();
    summary.metadata = metadata;
    summary.save_filename = save_filename;
    all_conditions_summary.(current_condition.name) = summary;
end % End of conditions loop

%% Final summary
total_elapsed = toc(overall_super_start_time);
fprintf('=== OVERALL PARAMETER SPACE EXPLORATION COMPLETE ===\n');
fprintf('Total duration: %.2f hours\n', total_elapsed/3600);
fprintf('Conditions analyzed: %d\n', length(conditions));
fprintf('==================================================\n');

% Save final summary
final_summary_filename = fullfile(output_dir_base, 'parameter_space_analysis_summary.mat');
save(final_summary_filename, 'all_conditions_summary', 'total_elapsed', 'conditions');
fprintf('Final summary of all conditions saved to %s\n', final_summary_filename);
fprintf('==================================================\n');

beep;pause(0.5);beep;pause(0.2);beep % done, wake up 

