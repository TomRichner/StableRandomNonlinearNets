% plot_parameter_space_results_paired_pulse_MI_multi.m
%
% Merge and plot consolidated parameter-space results for the paired-pulse
% MI study from MULTIPLE experiment folders. Each dataset can contribute
% only runs with LLE falling in a specified [L_min L_max] range.
%
% - Robust to missing conditions in any dataset
% - Aligns MI vs delay vectors across datasets by interpolation onto the
%   first encountered delay vector template
% - Produces the same set of plots as the single-dataset script

close all;
clear;
clc;

%% --- Script Configuration ---
counts_or_pdf = 'pdf'; % 'counts' or 'pdf'
invisible_y_axis = true;
use_bias_corrected_MI = true; % true for bias-corrected MI, false for uncorrected MI

% LLE histogram parameters
lle_range = [-1.5, 1.5];
n_bins_lle = 25;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins_lle), inf];

% LLE imagesc parameters
imagesc_lle_range = [-10 1]; % LLE range for the imagesc plot
n_bins_imagesc_lle_range = diff(imagesc_lle_range)+4;

% Mean rate histogram parameters
rate_range = [0, 5];
n_bins_rate = 25;
rate_bins = [linspace(rate_range(1), rate_range(2), n_bins_rate + 1), inf];

% MI parameters
mi_delay_for_histogram_samples = 200; % delay in samples for histogram column
delays_for_swarm_samples = [50]; % delays in samples for swarm plot
mi_range = [0, 3.5];
n_bins_mi = 25;
mi_bins = [linspace(mi_range(1), mi_range(2), n_bins_mi + 1), inf];

% LLE vs MI slice plot parameters
mi_delay_window_for_slice_samples = [150 250]; % delay window in samples

% Bootstrap parameters for confidence intervals
n_bootstrap_samples = 1000;
min_samples_for_ci = 5; % Minimum number of data points in a bin to compute CI

% Figure bin styling
outer_bin_width_multiplier = 1;
middle_ticks_lle = [0];
middle_ticks_rate = [0];
middle_ticks_mi = [0];

% Readable titles for conditions
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

%% Select multiple experiment folders and LLE ranges
prompt_num = inputdlg({'How many datasets to merge?'}, 'Multi-dataset setup', 1, {'2'});
if isempty(prompt_num)
    fprintf('User cancelled. Exiting.\n');
    return;
end
num_datasets = str2double(prompt_num{1});
if isnan(num_datasets) || num_datasets < 1
    fprintf('Invalid number. Exiting.\n');
    return;
end

datasets = struct('base_dir', {}, 'lle_range', {});
for di = 1:num_datasets
    fprintf('Select dataset folder %d of %d...\n', di, num_datasets);
    base_dir_i = uigetdir(pwd, sprintf('Select Experiment Directory (dataset %d)', di));
    if isequal(base_dir_i, 0)
        fprintf('User cancelled. Exiting.\n');
        return;
    end
    prompt_lr = inputdlg({sprintf('LLE range for dataset %d [L_min L_max] (e.g., 0 1):', di)}, ...
                         'LLE range', 1, {'0 1'});
    if isempty(prompt_lr)
        fprintf('User cancelled. Exiting.\n');
        return;
    end
    lr_vals = str2num(prompt_lr{1}); %#ok<ST2NM>
    if ~isnumeric(lr_vals) || numel(lr_vals) ~= 2
        fprintf('Invalid LLE range. Exiting.\n');
        return;
    end
    lr_vals = sort(lr_vals(:))';
    datasets(di).base_dir = base_dir_i;
    datasets(di).lle_range = lr_vals;
    fprintf('Dataset %d: %s with LLE range [%g, %g]\n', di, base_dir_i, lr_vals(1), lr_vals(2));
end

conditions_requested = {'no_adaptation', 'std_only', 'sfa_only', 'sfa_and_std'};

%% Aggregate data across datasets
all_data_combined = struct();
conditions_to_plot = {};

% Template for delay vector and fs
template_delay_vec = [];
template_fs = [];

% Helper function for MI alignment
align_mi_to_template = @(mi_vec, delay_vec, tmpl_delay) ...
    interp1(double(delay_vec(:)), double(mi_vec(:)), double(tmpl_delay(:)), 'linear', NaN)';

for di = 1:num_datasets
    base_dir = datasets(di).base_dir;
    lle_min = datasets(di).lle_range(1);
    lle_max = datasets(di).lle_range(2);
    fprintf('\nScanning dataset %d: %s\n', di, base_dir);

    for i = 1:length(conditions_requested)
        condition_name = conditions_requested{i};
        data_file = fullfile(base_dir, condition_name, sprintf('param_space_results_%s.mat', condition_name));
        if ~exist(data_file, 'file')
            fprintf('  - Missing condition "%s" in this dataset. Skipping.\n', condition_name);
            continue;
        end

        fprintf('  + Loading: %s\n', data_file);
        loaded = load(data_file);
        if ~isfield(loaded, 'results')
            fprintf('    ! No results field found. Skipping.\n');
            continue;
        end

        % Initialize combined storage for this condition
        if ~isfield(all_data_combined, condition_name)
            all_data_combined.(condition_name).lles = [];
            all_data_combined.(condition_name).mean_rates = [];
            all_data_combined.(condition_name).all_mi_vectors = {}; % rows interpolated/aligned to template
        end

        results = loaded.results;

        % Pass 1: Find a template delay vector if not set yet
        if isempty(template_delay_vec)
            for j = 1:numel(results)
                res = results{j};
                if isstruct(res) && isfield(res, 'success') && res.success && isfield(res, 'simulation_output')
                    so = res.simulation_output;
                    if use_bias_corrected_MI
                        mi_field_name = 'MI_vs_delay_corrected';
                    else
                        mi_field_name = 'MI_vs_delay_uncorrected';
                    end
                    if isfield(so, mi_field_name) && isfield(so, 'delay_vec')
                        delay_vec = so.delay_vec; mi_vec = so.(mi_field_name);
                        if isvector(delay_vec) && isvector(mi_vec) && numel(delay_vec) == numel(mi_vec)
                            template_delay_vec = delay_vec(:)';
                            if isfield(res, 'parameters') && isfield(res.parameters, 'fs')
                                template_fs = res.parameters.fs;
                            end
                            break;
                        end
                    end
                end
            end
        end

        % Pass 2: Collect and filter by dataset LLE range, align MI
        for j = 1:numel(results)
            res = results{j};
            if ~(isstruct(res) && isfield(res, 'success') && res.success)
                continue;
            end

            LLE_val = NaN; mean_rate_val = NaN;
            if isfield(res, 'simulation_output') && isfield(res.simulation_output, 'LLE')
                LLE_val = res.simulation_output.LLE;
            end
            if isfield(res, 'simulation_output') && isfield(res.simulation_output, 'mean_rate')
                mean_rate_val = res.simulation_output.mean_rate;
            end

            if isempty(LLE_val) || isnan(LLE_val) || LLE_val < lle_min || LLE_val > lle_max
                continue; % Outside desired LLE range for this dataset
            end

            % MI fields - choose between corrected and uncorrected
            if use_bias_corrected_MI
                mi_field_name = 'MI_vs_delay_corrected';
            else
                mi_field_name = 'MI_vs_delay_uncorrected';
            end
            so = res.simulation_output;
            if isfield(so, mi_field_name) && isfield(so, 'delay_vec')
                delay_vec = so.delay_vec; mi_vec = so.(mi_field_name);
                if ~isempty(template_delay_vec) && isvector(delay_vec) && isvector(mi_vec) && numel(delay_vec) == numel(mi_vec)
                    if isequal(delay_vec(:)', template_delay_vec(:)')
                        mi_vec_aligned = mi_vec(:)';
                    else
                        % Interpolate onto the template delay vector
                        try
                            mi_vec_aligned = align_mi_to_template(mi_vec, delay_vec, template_delay_vec);
                        catch
                            mi_vec_aligned = NaN(size(template_delay_vec));
                        end
                    end
                else
                    mi_vec_aligned = NaN(size(template_delay_vec));
                end
            else
                mi_vec_aligned = NaN(size(template_delay_vec));
            end

            % Append
            all_data_combined.(condition_name).lles(end+1) = LLE_val; %#ok<AGROW>
            all_data_combined.(condition_name).mean_rates(end+1) = mean_rate_val; %#ok<AGROW>
            all_data_combined.(condition_name).all_mi_vectors{end+1} = mi_vec_aligned; %#ok<AGROW>
            if ~ismember(condition_name, conditions_to_plot)
                conditions_to_plot{end+1} = condition_name; %#ok<AGROW>
            end
        end
    end
end

% Validate combined data
num_conditions = length(conditions_to_plot);
if num_conditions == 0
    fprintf('No valid condition data found across selected datasets. Exiting.\n');
    return;
end

fprintf('\nDatasets merged. Found conditions: %s\n', strjoin(conditions_to_plot, ', '));
if isempty(template_delay_vec)
    fprintf('Warning: No valid MI delay vectors found. MI plots involving delay will be empty.\n');
end

%% Extract values (already aggregated)
extracted_values = struct();
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    extracted_values.(condition_name).lles = all_data_combined.(condition_name).lles;
    extracted_values.(condition_name).mean_rates = all_data_combined.(condition_name).mean_rates;
    extracted_values.(condition_name).all_mi_vectors = all_data_combined.(condition_name).all_mi_vectors;
    extracted_values.(condition_name).delay_vec = template_delay_vec;
    extracted_values.(condition_name).fs = template_fs;
    % For MI hist and swarm: derive from aligned vectors if available
    mi_for_hist = [];
    mi_for_swarm = [];
    mi_vectors = extracted_values.(condition_name).all_mi_vectors;
    if ~isempty(template_delay_vec)
        [~, delay_idx_hist] = min(abs(template_delay_vec - mi_delay_for_histogram_samples));
        delay_indices_swarm = zeros(1, length(delays_for_swarm_samples));
        for k_delay = 1:length(delays_for_swarm_samples)
            [~, delay_indices_swarm(k_delay)] = min(abs(template_delay_vec - delays_for_swarm_samples(k_delay)));
        end
        for k = 1:length(mi_vectors)
            v = mi_vectors{k};
            if ~isempty(v) && numel(v) == numel(template_delay_vec)
                if ~isnan(v(delay_idx_hist))
                    mi_for_hist(end+1) = v(delay_idx_hist); %#ok<AGROW>
                end
                vals = v(delay_indices_swarm);
                mi_for_swarm(end+1, :) = vals; %#ok<AGROW>
            end
        end
    end
    extracted_values.(condition_name).mi_for_hist = mi_for_hist;
    extracted_values.(condition_name).mi_for_swarm = mi_for_swarm;
    fprintf('Condition "%s": %d LLEs, %d mean rates, %d full MI vectors.\n', ...
        condition_name, length(extracted_values.(condition_name).lles), ...
        length(extracted_values.(condition_name).mean_rates), ...
        length(extracted_values.(condition_name).all_mi_vectors));
end

%% Create main distributions figure (LLE and rate)
if use_bias_corrected_MI
    mi_type_str = ' (Bias-Corrected MI)';
else
    mi_type_str = ' (Uncorrected MI)';
end
fig_main = figure('Name', ['LLE and Rate Distributions (Paired-Pulse, Multi)' mi_type_str], 'Position', [100, 100, 650, 1020]);

% Column 1: LLE
ax_lle = gobjects(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_lle(i) = subplot(num_conditions, 2, 2*i-1);
    [counts, ~] = histcounts(extracted_values.(condition_name).lles, lle_bins);
    % Finite bins for LLE (two-sided outer bins)
    std_width_lle = (lle_range(2) - lle_range(1)) / (n_bins_lle - 1);
    middle_edges_lle = linspace(lle_range(1), lle_range(2), n_bins_lle);
    first_edge = middle_edges_lle(1) - outer_bin_width_multiplier * std_width_lle;
    last_edge = middle_edges_lle(end) + outer_bin_width_multiplier * std_width_lle;
    lle_bins_finite = [first_edge, middle_edges_lle, last_edge];
    if strcmpi(counts_or_pdf, 'pdf')
        plot_values = counts ./ (sum(counts) + eps);
    else
        plot_values = counts;
    end
    histogram('BinEdges', lle_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    if i == num_conditions, xlabel('LLE (\lambda_1)','FontSize',22); end
    if invisible_y_axis, set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none'); end
    box off;
    ax = gca; xlim([lle_bins_finite(1), lle_bins_finite(end)]);
    % ticks
    center_first_bin = (lle_bins_finite(1) + lle_bins_finite(2)) / 2;
    center_last_bin  = (lle_bins_finite(end-1) + lle_bins_finite(end)) / 2;
    final_ticks = [center_first_bin, middle_ticks_lle, center_last_bin];
    new_labels = [sprintf('<= %.1f', lle_range(1)), arrayfun(@(x) sprintf('%.1f', x), middle_ticks_lle, 'UniformOutput', false), {sprintf('>= %.1f', lle_range(2))}];
    xticks(ax, final_ticks); xticklabels(ax, new_labels); ax.XTickLabelRotation = 45;
    % ylabel as text on each row
    if isKey(condition_titles, condition_name)
        title_str = condition_titles(condition_name);
    else
        title_str = condition_name;
    end
    text(-0.2, 0.5, sprintf('%s', title_str), 'Units','normalized', 'Rotation',90, 'FontSize',22, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

% Column 2: Mean Rate
ax_rate = gobjects(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_rate(i) = subplot(num_conditions, 2, 2*i);
    [counts, ~] = histcounts(extracted_values.(condition_name).mean_rates, rate_bins);
    % Finite bins for Rate (lower bound at zero)
    std_width_rate = (rate_range(2) - rate_range(1)) / n_bins_rate;
    middle_edges_rate = linspace(rate_range(1), rate_range(2), n_bins_rate + 1);
    first_edge = middle_edges_rate(1);
    last_edge = middle_edges_rate(end) + outer_bin_width_multiplier * std_width_rate;
    rate_bins_finite = [first_edge, middle_edges_rate(2:end), last_edge];
    if strcmpi(counts_or_pdf, 'pdf')
        plot_values = counts ./ (sum(counts) + eps);
    else
        plot_values = counts;
    end
    histogram('BinEdges', rate_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    if i == num_conditions, xlabel('Mean Firing Rate (Hz)','FontSize',22); end
    if invisible_y_axis, set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none'); end
    box off;
    ax = gca; xlim([rate_bins_finite(1), rate_bins_finite(end)]);
    center_last_bin = (rate_bins_finite(end-1) + rate_bins_finite(end)) / 2;
    final_ticks = [middle_ticks_rate, center_last_bin];
    new_labels = [arrayfun(@(x) sprintf('%.0f', x), middle_ticks_rate, 'UniformOutput', false), {sprintf('>= %.0f', rate_range(2))}];
    xticks(ax, final_ticks); xticklabels(ax, new_labels); ax.XTickLabelRotation = 45;
end

% Link y axes
linkaxes([ax_lle; ax_rate], 'y');

%% Imagesc figure: Collapse across conditions - LLE vs delay with MI
pooled_delay = [];
pooled_mi = [];
pooled_lle = [];

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    delay_vec = extracted_values.(condition_name).delay_vec;
    mi_vectors = extracted_values.(condition_name).all_mi_vectors;
    lle_list = extracted_values.(condition_name).lles;
    if isempty(delay_vec) || isempty(mi_vectors)
        continue
    end
    for k = 1:length(mi_vectors)
        mi_vec = mi_vectors{k};
        if isempty(mi_vec) || ~isvector(mi_vec)
            continue
        end
        if length(mi_vec) ~= length(delay_vec)
            continue
        end
        pooled_mi(end+1, :) = mi_vec(:)'; %#ok<AGROW>
        if ~isempty(lle_list)
            lle_val = lle_list(min(k, length(lle_list)));
        else
            lle_val = NaN;
        end
        pooled_lle(end+1, 1) = lle_val; %#ok<AGROW>
    end
end

fig_img = figure('Name', ['MI vs Delay by LLE (Collapsed Across Conditions, Multi)' mi_type_str], 'Position', [643  500  404   364]);
if isempty(pooled_mi) || isempty(template_delay_vec)
    axes; axis off; text(0.5,0.5,'No MI data available for imagesc','HorizontalAlignment','center');
else
    % Bin LLE into edges and compute median MI in each LLE bin per delay
    lle_edges = linspace(imagesc_lle_range(1), imagesc_lle_range(2), n_bins_imagesc_lle_range);
    n_lle_bins = length(lle_edges)-1;
    n_delays = length(template_delay_vec);
    mi_img = NaN(n_lle_bins, n_delays);

    for b = 1:n_lle_bins
        in_bin = pooled_lle >= lle_edges(b) & pooled_lle < lle_edges(b+1);
        if any(in_bin)
            mi_in_bin = squeeze(pooled_mi(in_bin, :));
            mi_img(b, :) = mean(mi_in_bin, 1, 'omitnan');
        else
            mi_img(b, :) = nan(1, size(pooled_mi, 2));
        end
    end

    imagesc((template_delay_vec./template_fs), (lle_edges(1:end-1)+lle_edges(2:end))/2, mi_img);
    axis xy; colormap(hot);
    xlabel('Delay (s)','FontSize',22);
    ylabel('LLE (binned center)','FontSize',22);
    title('Mean paired-pulse MI vs delay and LLE (multi)')
    colorbar; ylabel(colorbar, 'MI (bits)');
end

%% LLE vs MI slice figure
fig_lle_mi_slice = figure('Name', ['LLE vs MI Slice (Multi)' mi_type_str], 'Position', [700, 400, 380 300]);
if isempty(pooled_mi) || isempty(template_delay_vec)
    axes; axis off; text(0.5,0.5,'No MI data available for LLE vs MI slice plot','HorizontalAlignment','center');
else
    % Find delay indices for the averaging window
    delay_indices = find(template_delay_vec >= mi_delay_window_for_slice_samples(1) & ...
                         template_delay_vec <= mi_delay_window_for_slice_samples(2));

    if isempty(delay_indices)
        axes; axis off; text(0.5,0.5,'Specified delay window is empty or out of bounds','HorizontalAlignment','center');
    else
        % Average MI over the specified delay window
        mi_averaged = mean(pooled_mi(:, delay_indices), 2, 'omitnan');
        
        % Remove runs where LLE or averaged MI is NaN
        valid_indices = ~isnan(pooled_lle) & ~isnan(mi_averaged);
        lle_valid = pooled_lle(valid_indices);
        
        % Optional: histogram of valid LLEs (separate figure)
        fig_lle_hist = figure('Name', 'LLE Distribution (Multi)', 'Position', [800, 300, 400, 300]); %#ok<NASGU>
        histogram(lle_valid, 100);
        xlabel('LLE'); ylabel('Count'); title('Distribution of Valid LLE Values (Multi)'); box off;
        figure(fig_lle_mi_slice);

        mi_valid = mi_averaged(valid_indices);

        % Use the same LLE bins as the imagesc plot for consistency
        lle_edges_slice = linspace(imagesc_lle_range(1), imagesc_lle_range(2), n_bins_imagesc_lle_range);
        n_lle_bins_slice = length(lle_edges_slice) - 1;
        
        lle_bin_centers = (lle_edges_slice(1:end-1) + lle_edges_slice(2:end)) / 2;
        mi_mean_in_bin = NaN(1, n_lle_bins_slice);
        mi_ci_lower = NaN(1, n_lle_bins_slice);
        mi_ci_upper = NaN(1, n_lle_bins_slice);
        
        for b = 1:n_lle_bins_slice
            in_bin_indices = find(lle_valid >= lle_edges_slice(b) & lle_valid < lle_edges_slice(b+1));
            if ~isempty(in_bin_indices)
                mi_in_bin = mi_valid(in_bin_indices);
                mi_mean_in_bin(b) = mean(mi_in_bin, 'omitnan');
                if length(mi_in_bin) >= min_samples_for_ci
                    ci = bootci(n_bootstrap_samples, @mean, mi_in_bin);
                    mi_ci_lower(b) = ci(1);
                    mi_ci_upper(b) = ci(2);
                end
            end
        end
        
        % Plotting with asymmetric confidence intervals
        y_neg = mi_mean_in_bin - mi_ci_lower;
        y_pos = mi_ci_upper - mi_mean_in_bin;
        errorbar(lle_bin_centers, mi_mean_in_bin, y_neg, y_pos, 'ko-', 'LineWidth', 1.5, 'CapSize', 4);
        xlabel('LLE'); ylabel('Mutual Information (bits)'); box off;
        set(gca, 'YTick', [0 1 2 3], 'XTick', -15:5:0);
    end
end

%% Create Swarm Plot figure for MI at different delays
fig_swarm = figure('Name', ['MI Swarm Plot by Condition and Delay (Multi)' mi_type_str], 'Position', [200, 200, 1200, 600]);
tiledlayout(1, num_conditions, 'TileSpacing', 'compact');
all_ax_swarm = [];
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    mi_data = extracted_values.(condition_name).mi_for_swarm;

    if isempty(mi_data)
        nexttile; axis off; title_str = condition_name; if isKey(condition_titles, condition_name), title_str = condition_titles(condition_name); end
        title(title_str); text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        continue;
    end

    ax_swarm = nexttile; all_ax_swarm = [all_ax_swarm, ax_swarm]; %#ok<AGROW>
    num_points = size(mi_data, 1); %#ok<NASGU>
    num_delays = length(delays_for_swarm_samples); %#ok<NASGU>
    mi_flat = mi_data(:);
    delay_labels_categorical = categorical(repmat(delays_for_swarm_samples', size(mi_data,1), 1));
    swarmchart(ax_swarm, delay_labels_categorical, mi_flat, 'filled');
    title_str = condition_name; if isKey(condition_titles, condition_name), title_str = condition_titles(condition_name); end
    title(title_str);
    xlabel('Delay (samples)');
    if i == 1
        ylabel('Mutual Information (bits)');
    else
        yticklabels({});
    end
end
if ~isempty(all_ax_swarm)
    linkaxes(all_ax_swarm, 'y');
end

%% Create Violin Plot figure for MI at different delays (using violinPlots2)
fig_violin = figure('Name', ['MI Violin Plot by Condition and Delay (Multi)' mi_type_str], 'Position', [200   354   610  364]);
tiledlayout(1, num_conditions, 'TileSpacing', 'compact');
all_ax_violin = [];

% Add the violinplot path if it's not already there
violinplot_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'src', 'supporting_functions', 'external', 'violinPlots2');
if ~exist('Violin.m', 'file')
    fprintf('Adding violinplot path to MATLAB path: %s\n', violinplot_path);
    addpath(violinplot_path);
end

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    mi_data = extracted_values.(condition_name).mi_for_swarm;

    ax_violin = nexttile; all_ax_violin = [all_ax_violin, ax_violin]; %#ok<AGROW>
    if isempty(mi_data)
        axis off; title_str = condition_name; if isKey(condition_titles, condition_name), title_str = condition_titles(condition_name); end
        title(title_str); text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        continue;
    end
    
    % Calculate bandwidth for this condition (20% of range)
    data_range = range(mi_data(:));
    bandwidth = 0.2 * data_range;
    
    delay_labels = string(delays_for_swarm_samples);
    violins = violinplot(mi_data, delay_labels, 'Parent', ax_violin, ...
        'Bandwidth', bandwidth, ...
        'KSDensityOptions', {'support', [-eps, 5]}, ...
        'ShowMedian', false, ...
        'QuartileStyle', 'shadow', ...
        'ShowBox', false, ...
        'ShowWhiskers', false, ...
        'MarkerSize', 9);
    hold(ax_violin, 'on');

    % Medians across delays
    medians = median(mi_data, 1, 'omitnan');
    for k = 1:length(violins)
        v = violins(k);
        y_median = medians(k);
        if isempty(v.ViolinPlot) || ~isvalid(v.ViolinPlot)
            continue;
        end
        violin_patch = v.ViolinPlot; x_data = violin_patch.XData; y_data = violin_patch.YData;
        [y_max, max_idx] = max(y_data);
        y_up = y_data(1:max_idx); x_up = x_data(1:max_idx);
        y_down = y_data(max_idx:end); x_down = x_data(max_idx:end);
        [y_up_unique, ia_up] = unique(y_up); x_up_unique = x_up(ia_up);
        [y_down_unique, ia_down] = unique(y_down); x_down_unique = x_down(ia_down);
        x1 = interp1(y_up_unique, x_up_unique, y_median, 'linear');
        x2 = interp1(y_down_unique, x_down_unique, y_median, 'linear');
        if isfinite(x1) && isfinite(x2)
            plot(ax_violin, [x1, x2], [y_median, y_median], 'k-', 'LineWidth', 2);
        end
    end
    hold(ax_violin, 'off');
    box(ax_violin, 'off');
    title_str = condition_name; if isKey(condition_titles, condition_name), title_str = condition_titles(condition_name); end
    xlabel(ax_violin, title_str);
    set(ax_violin, 'XTick', []);
    if i == 1
        ylabel('Mutual Information (bits)'); yticks(ax_violin, [0 1 2 3]);
    else
        set(ax_violin, 'YTick', [], 'YColor', 'none'); ylabel('');
    end
end
if ~isempty(all_ax_violin)
    linkaxes(all_ax_violin, 'y'); ylim(all_ax_violin, [0 2.55]);
end

%% Save figures
output_dir_for_figs = fullfile(pwd, 'analysis_plots_multi');
if ~exist(output_dir_for_figs, 'dir')
    mkdir(output_dir_for_figs);
end
if use_bias_corrected_MI
    mi_suffix = '_bias_corrected';
else
    mi_suffix = '_uncorrected';
end
save_some_figs_to_folder_2(output_dir_for_figs, ['PPMI_MULTI_LLE_rate_MI_distributions_v1' mi_suffix], fig_main.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, ['PPMI_MULTI_MI_vs_delay_LLE_imagesc_v1' mi_suffix], fig_img.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, ['PPMI_MULTI_LLE_vs_MI_slice_v1' mi_suffix], fig_lle_mi_slice.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, ['PPMI_MULTI_MI_swarm_by_condition_v1' mi_suffix], fig_swarm.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, ['PPMI_MULTI_MI_violin_by_condition_v1' mi_suffix], fig_violin.Number, {'fig', 'svg', 'png'});

fprintf('Multi-dataset plotting complete. Figures saved to %s\n', output_dir_for_figs);


