% scatter_LLE_vs_rate_jitter_mixed.m
%
% This script loads the fully consolidated data from a parameter space search
% and creates a single scatter plot of LLE vs. mean spike rate, with
% jitter added to each point to prevent overplotting.
% The data from all conditions is mixed and plotted together to avoid visual
% overrepresentation of any single condition.
% Each point is a network simulation, color-coded by its adaptation condition.

close all;
clear all;
clc;

%% --- Script Configuration ---

% Readable titles and colors for conditions
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

condition_colors = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]}); % Default MATLAB colors

% Max values for filtering
max_lle = 5;
max_rate = 10;
min_lle = -10;

%% Get the base directory for the analysis
fprintf('Please select the main experiment directory to plot...\n');
base_dir = uigetdir(pwd, 'Select Experiment Directory');
if isequal(base_dir, 0)
    fprintf('User cancelled. Exiting.\n');
    return;
end
fprintf('Selected directory: %s\n', base_dir);

%% Find and load the data for each condition
conditions_to_plot = {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'};
num_conditions = length(conditions_to_plot);
all_data = struct();
found_all = true;

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    data_file = fullfile(base_dir, condition_name, sprintf('param_space_results_%s.mat', condition_name));
    
    if exist(data_file, 'file')
        fprintf('Loading: %s\n', data_file);
        all_data.(condition_name) = load(data_file);
    else
        fprintf('Error: Data file not found for condition "%s".\nExpected at: %s\n', condition_name, data_file);
        found_all = false;
    end
end

if ~found_all
    fprintf('Could not find all required data files. Exiting.\n');
    return;
end

%% Process the data to extract LLE and mean rates for plotting
extracted_values = struct();
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    results = all_data.(condition_name).results;
    
    lles = [];
    mean_rates = [];
    
    num_results = length(results);
    for j = 1:num_results
        res = results{j};
        if isstruct(res) && isfield(res, 'success') && res.success && ...
           isfield(res.simulation_output, 'LLE') && isfield(res.simulation_output, 'mean_rate')
            
            lle_val = res.simulation_output.LLE;
            rate_val = res.simulation_output.mean_rate;
            
            % Ensure values are single scalars and not empty or NaN, then clamp
            if isscalar(lle_val) && ~isnan(lle_val) && isscalar(rate_val) && ~isnan(rate_val)
                
                % Clamp LLE to be within [min_lle, max_lle]
                lle_val = max(min(lle_val, max_lle), min_lle);
                
                % Clamp rate to be no more than max_rate
                rate_val = min(rate_val, max_rate);

                 lles(end+1) = lle_val;
                 mean_rates(end+1) = rate_val;
            end
        end
    end
    extracted_values.(condition_name).lles = lles;
    extracted_values.(condition_name).mean_rates = mean_rates;
    fprintf('Condition "%s": Found %d valid (LLE, rate) pairs.\n', ...
        condition_name, length(lles));
end

%% --- Combine and Shuffle Data for Plotting ---
all_lles = [];
all_mean_rates = [];
all_colors = [];

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    
    lles = extracted_values.(condition_name).lles;
    mean_rates = extracted_values.(condition_name).mean_rates;
    color = condition_colors(condition_name);
    num_points = length(lles);
    
    all_lles = [all_lles, lles];
    all_mean_rates = [all_mean_rates, mean_rates];
    all_colors = [all_colors; repmat(color, num_points, 1)];
end

% Shuffle the data to mix conditions during plotting
num_total_points = length(all_lles);
shuffle_idx = randperm(num_total_points);

all_lles_shuffled = all_lles(shuffle_idx);
all_mean_rates_shuffled = all_mean_rates(shuffle_idx);
all_colors_shuffled = all_colors(shuffle_idx, :);


%% Create the scatter plot
fig_scatter = figure('Name', 'LLE vs. Mean Rate', 'Position', [100, 100, 584,  446]);
hold on;

% Add jitter to all data points
jitter_std = 0.04;
lles_jittered = all_lles_shuffled + jitter_std * randn(size(all_lles_shuffled));
mean_rates_jittered = all_mean_rates_shuffled + jitter_std * randn(size(all_mean_rates_shuffled));

% Plot all points in one go with specified colors
scatter(lles_jittered, mean_rates_jittered, 14, all_colors_shuffled, 'filled', 'MarkerFaceAlpha', 0.4);

% Create dummy plots for the legend
plot_handles = gobjects(num_conditions, 1);
condition_labels = cell(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    color = condition_colors(condition_name);
    title_str = condition_titles(condition_name);
    
    plot_handles(i) = scatter(NaN, NaN, 50, color, 'filled');
    condition_labels{i} = title_str;
end

hold off;

% Add plot labels and title
xlabel('$\lambda_1$', 'Interpreter', 'latex','FontSize',28);
ylabel('Mean Firing Rate (Hz)');

% Custom Rate Ticks
min_rate = 0; % Assuming minimum rate is 0
middle_ticks_rate = [5];
final_ticks_rate = [min_rate, middle_ticks_rate, max_rate];

% Create labels for the ticks
new_labels_rate = cell(1, length(final_ticks_rate));
new_labels_rate{1} = sprintf('%.0f', min_rate);
for k = 1:length(middle_ticks_rate)
    new_labels_rate{1+k} = sprintf('%.0f', middle_ticks_rate(k));
end
new_labels_rate{end} = sprintf('>= %.0f', max_rate);

ax = gca;
yticks(ax, final_ticks_rate);
yticklabels(ax, new_labels_rate);

% Custom LLE Ticks
middle_ticks_lle = [-5, 0];
final_ticks = [min_lle, middle_ticks_lle, max_lle];

% Create labels for the ticks
new_labels = cell(1, length(final_ticks));
new_labels{1} = sprintf('<= %.1f', min_lle);
for k = 1:length(middle_ticks_lle)
    new_labels{1+k} = sprintf('%.1f', middle_ticks_lle(k));
end
new_labels{end} = sprintf('>= %.1f', max_lle);

ax = gca;
xticks(ax, final_ticks);
xticklabels(ax, new_labels);
ax.XTickLabelRotation = 45;

% Add a legend
legend(plot_handles, condition_labels, 'Location', 'northeastoutside', 'FontSize', 14);

% Improve axis appearance

% Add a vertical line at LLE = 0
h = xline(0, '--k', 'LineWidth', 1.5);
h.HandleVisibility = 'off';
xlim([min_lle max_lle]+[-0.1 0.1])
ylim([0 max_rate]+[-0.1 0.15])

%% Save the figure
output_dir_for_figs = fullfile(base_dir, 'analysis_plots_jitter_mixed');

fprintf('\nSaving figure to: %s\n', output_dir_for_figs);
% Assuming save_some_figs_to_folder_2 is on the path
% If not, you might need to add its location to the path.
save_some_figs_to_folder_2(output_dir_for_figs, 'LLE_vs_rate_scatter_jitter_mixed', fig_scatter.Number, {'fig', 'svg', 'png'});

fprintf('Plotting complete.\n');
beep; pause(0.5); beep % wake up
