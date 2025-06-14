% Script for sensitivity analysis of the SRNN model

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



params_to_analyze = {'EE_factor','mean_weight','n'}; % Only analyze these parameters
% params_to_analyze = {}; % Empty means analyze all parameters

%% Analysis Parameters
n_levels = 5; % Number of values to test for each parameter
n_reps = 5;   % Number of repetitions with different random seeds for each level

note = 'quick_LLE_and_SR'

% Timestamp for folder name
dt_str = lower(strrep(datestr(now, 'mmm_dd_yy_hh_MM_AM'), ':', '_'));

%% Create an ABSOLUTE base directory and be sure it exists
output_dir_base = fullfile(pwd, ...
    ['sensi_' note '_nLevs_' num2str(n_levels) '_nReps_' num2str(n_reps) '_' dt_str]);

if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);      % create the directory tree once, up-front
end

%% after the output_dir_base is made, copy current mfile (this file) into it
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
    
    %% Default Simulation Parameters (based on SRNN_caller.m)
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
    p_default.tau_b_E_2 = 1; % if n_b_E == 1, then this value is used.
    p_default.tau_STD = 0.5;
    p_default.c_SFA_factor = 0.5;
    p_default.n_a_E = current_condition.n_a_E_val;
    p_default.n_b_E = current_condition.n_b_E_val;
    
    %% Parameter Ranges for Sensitivity Analysis
    ranges.fs = [250, 2000];
    ranges.n = [10, 50];
    ranges.EE_factor = [0, 4];
    ranges.IE_factor = [0, 4];
    ranges.EI = [1/p_default.n, 1.0];  % Changed from [0, 1.0] to avoid EI=0 issues
    ranges.E_self = [0.0, 0.5];
    ranges.mean_weight = [1/p_default.n, 2];
    ranges.DC = [0, 4];
    ranges.mean_in_out_degree = [1, p_default.n-1];
    ranges.tau_a_E_2 = [1, 30];
    ranges.tau_b_E_2 = [2, 30];
    ranges.tau_STD = [0.2, 1];
    ranges.c_SFA_factor = [0.0, 4.0];
    ranges.n_a_E = [0, 10];
    ranges.n_b_E = [0, 5];
    
    param_names = fieldnames(ranges);
    
    % Exclude condition parameters from sensitivity analysis
    param_names = setdiff(param_names, {'n_a_E', 'n_b_E'}, 'stable');
    
    % Filter to only analyze specified parameters
    if ~isempty(params_to_analyze)
        % Validate that all specified parameters exist in ranges
        invalid_params = setdiff(params_to_analyze, param_names);
        if ~isempty(invalid_params)
            error('Invalid parameter names specified: %s', strjoin(invalid_params, ', '));
        end
        
        % Filter param_names to only include specified parameters
        param_names = intersect(param_names, params_to_analyze, 'stable');
        fprintf('Running sensitivity analysis for subset of parameters: %s\n', strjoin(param_names, ', '));
    else
        fprintf('Running sensitivity analysis for all parameters\n');
    end
    
    output_dir = fullfile(output_dir_base, current_condition.name);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Initialize overall timing and resource monitoring
    overall_start_time = tic;
    total_params = length(param_names);
    total_runs_all = total_params * n_levels * n_reps;
    
    fprintf('=== SRNN Sensitivity Analysis Started for %s ===\n', current_condition.name);
    fprintf('Total parameters: %d\n', total_params);
    fprintf('Levels per parameter: %d\n', n_levels);
    fprintf('Repetitions per level: %d\n', n_reps);
    fprintf('Total runs: %d\n', total_runs_all);
    fprintf('==========================================\n\n');

    %% Run Sensitivity Analysis for each parameter
    all_results_summary = struct();

    for i = 1:length(param_names)
        param_name = param_names{i};
        param_range = ranges.(param_name);

        fprintf('--- Starting sensitivity analysis for: %s (%d/%d) ---\n', param_name, i, total_params);
        
        % Resource monitoring
        try
            if ispc
                [~, mem_info] = system('wmic OS get TotalVisibleMemorySize,FreePhysicalMemory /value');
                fprintf('System memory status: %s\n', strtrim(mem_info));
            elseif ismac || isunix
                [~, mem_info] = system('vm_stat | grep "Pages free\|Pages active\|Pages inactive"');
                fprintf('Memory status:\n%s\n', mem_info);
            end
        catch
            fprintf('Could not retrieve memory status\n');
        end
        
        param_start_time = tic;

        % Create vector of parameter levels
        if ismember(param_name, {'n', 'n_a_E', 'mean_in_out_degree'})
            param_levels = round(linspace(param_range(1), param_range(2), n_levels));
        else
            param_levels = linspace(param_range(1), param_range(2), n_levels);
        end

        total_runs = n_levels * n_reps;
        results = cell(total_runs, 1);
        
        % Track success/failure statistics
        n_success = 0;
        n_failed = 0;
        error_types = {};

        parfor k = 1:total_runs
        % for k = 1:total_runs
            level_idx = floor((k-1) / n_reps) + 1;
            rep_idx = mod(k-1, n_reps) + 1;
            
            current_p = p_default;
            current_p.(param_name) = param_levels(level_idx);
            
            sim_seed = k + (i-1)*total_runs; % Ensure unique seeds across all parameters
            
            % Enhanced progress reporting
            if mod(k, max(1, floor(total_runs/20))) == 1 || k <= 5 || k == total_runs
                fprintf('Running %s level %d/%d (val: %.3f), rep %d/%d (run %d/%d) [seed=%d]\n', ...
                    param_name, level_idx, n_levels, param_levels(level_idx), rep_idx, n_reps, k, total_runs, sim_seed);
            end

            run_start = tic;
            try
                result = SRNN_caller_wrapped_for_sensitivity_dual_stage(...
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
                result.run_duration = run_duration;
                result.success = true;
                result.param_value = param_levels(level_idx);
                result.seed = sim_seed;
                results{k} = result;
                
            catch ME
                run_duration = toc(run_start);
                fprintf('ERROR in %s level %d (val: %.3f), rep %d [run %d/%d, seed=%d]:\n', ...
                    param_name, level_idx, param_levels(level_idx), rep_idx, k, total_runs, sim_seed);
                fprintf('  Error: %s\n', ME.message);
                if ~isempty(ME.stack)
                    fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
                end
            
                % Store detailed error information
                results{k} = struct(...
                    'success', false, ...
                    'error', ME, ...
                    'error_message', ME.message, ...
                    'param_value', param_levels(level_idx), ...
                    'seed', sim_seed, ...
                    'run_duration', run_duration, ...
                    'level_idx', level_idx, ...
                    'rep_idx', rep_idx ...
                );
            end
        end
        
        % Calculate statistics
        success_mask = cellfun(@(x) isfield(x, 'success') && x.success, results);
        n_success = sum(success_mask);
        n_failed = total_runs - n_success;
        
        if n_failed > 0
            failed_results = results(~success_mask);
            error_messages = cellfun(@(x) x.error_message, failed_results, 'UniformOutput', false);
            [unique_errors, ~, error_idx] = unique(error_messages);
            error_counts = accumarray(error_idx, 1);
            
            fprintf('\n--- Error Summary for %s ---\n', param_name);
            for j = 1:length(unique_errors)
                fprintf('  "%s": %d occurrences\n', unique_errors{j}, error_counts(j));
            end
        end
        
        elapsed_time = toc(param_start_time);
        fprintf('--- Analysis for %s finished in %.2f minutes ---\n', param_name, elapsed_time/60);
        fprintf('Success rate: %d/%d (%.1f%%)\n', n_success, total_runs, 100*n_success/total_runs);

        % Reshape results array to be n_levels x n_reps
        results_reshaped = reshape(results, [n_reps, n_levels])';

        % Save results with enhanced metadata
        save_filename = fullfile(output_dir, sprintf('sensitivity_%s.mat', param_name));
        metadata = struct();
        metadata.param_name = param_name;
        metadata.param_levels = param_levels;
        metadata.param_range = param_range;
        metadata.p_default = p_default;
        metadata.n_reps = n_reps;
        metadata.n_levels = n_levels;
        metadata.elapsed_time = elapsed_time;
        metadata.n_success = n_success;
        metadata.n_failed = n_failed;
        metadata.success_rate = n_success/total_runs;
        metadata.analysis_date = datestr(now);
        metadata.condition = current_condition;
        
        save(save_filename, 'results_reshaped', 'metadata' ,'-v7.3');
        fprintf('Results saved to %s\n', save_filename);
        
        % Store summary for overall analysis
        all_results_summary.(param_name) = metadata;
        
        % Intermediate backup every 3 parameters
        if mod(i, 3) == 0
            backup_filename = fullfile(output_dir, sprintf('intermediate_backup_params_%d.mat', i));
            save(backup_filename, 'all_results_summary');
            fprintf('Intermediate backup saved to %s\n', backup_filename);
        end
        
        % Clear large variables to free memory
        clear results results_reshaped;
        
        % Estimate remaining time
        if i < total_params
            avg_time_per_param = elapsed_time;
            if i > 1
                % Use average of completed parameters
                completed_times = zeros(1, i);
                for j = 1:i
                    p_name = param_names{j};
                    completed_times(j) = all_results_summary.(p_name).elapsed_time;
                end
                avg_time_per_param = mean(completed_times);
            end
            remaining_time_est = avg_time_per_param * (total_params - i);
            fprintf('Estimated remaining time: %.1f hours\n', remaining_time_est/3600);
        end
        
        fprintf('\n');
    end

    all_conditions_summary.(current_condition.name) = all_results_summary;
end % End of conditions loop

%% Final summary
total_elapsed = toc(overall_super_start_time);
fprintf('=== OVERALL SENSITIVITY ANALYSIS COMPLETE ===\n');
fprintf('Total duration: %.2f hours\n', total_elapsed/3600);
fprintf('Conditions analyzed: %d\n', length(conditions));

% Calculate overall statistics
total_success_all_conditions = 0;
total_runs_all_conditions = 0;

condition_names_processed = fieldnames(all_conditions_summary);
for c_idx = 1:length(condition_names_processed)
    condition_name = condition_names_processed{c_idx};
    summary_for_condition = all_conditions_summary.(condition_name);
    param_names_in_summary = fieldnames(summary_for_condition);
    
    total_success_condition = 0;
    total_runs_condition = 0;
    
    for p_idx = 1:length(param_names_in_summary)
        param_name = param_names_in_summary{p_idx};
        if isfield(summary_for_condition, param_name)
            total_success_condition = total_success_condition + summary_for_condition.(param_name).n_success;
            total_runs_condition = total_runs_condition + (summary_for_condition.(param_name).n_success + summary_for_condition.(param_name).n_failed);
        end
    end
    
    if total_runs_condition > 0
        fprintf('  Condition: %s, Success rate: %d/%d (%.1f%%)\n', ...
            condition_name, total_success_condition, total_runs_condition, 100*total_success_condition/total_runs_condition);
    else
        fprintf('  Condition: %s, No runs completed.\n', condition_name);
    end
        
    total_success_all_conditions = total_success_all_conditions + total_success_condition;
    total_runs_all_conditions = total_runs_all_conditions + total_runs_condition;
end

if total_runs_all_conditions > 0
    overall_success_rate_val = 100*total_success_all_conditions/total_runs_all_conditions;
    fprintf('Overall success rate across all conditions: %d/%d (%.1f%%)\n', ...
        total_success_all_conditions, total_runs_all_conditions, overall_success_rate_val);
else
    overall_success_rate_val = 0;
    fprintf('No runs completed across all conditions.\n');
end

% Save final summary
final_summary_filename = fullfile(output_dir_base, 'sensitivity_analysis_summary_all_conditions.mat');
summary_data = struct();
summary_data.all_conditions_summary = all_conditions_summary;
summary_data.total_elapsed_time = total_elapsed;
summary_data.analysis_completed = datestr(now);
if total_runs_all_conditions > 0
    summary_data.overall_success_rate = total_success_all_conditions/total_runs_all_conditions;
else
    summary_data.overall_success_rate = 0;
end
summary_data.conditions = conditions;

save(final_summary_filename, 'summary_data', '-v7.3');
fprintf('Final summary saved to %s\n', final_summary_filename);
fprintf('======================================\n');

beep;pause(0.5);beep;pause(0.2);beep % done, wake up