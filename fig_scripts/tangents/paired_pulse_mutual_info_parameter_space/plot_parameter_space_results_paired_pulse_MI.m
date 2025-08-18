% plot_parameter_space_results_paired_pulse_MI.m
%
% Loads consolidated parameter-space results for the paired-pulse MI study
% and creates distributions of LLE, rate, MI at a chosen delay, as well as
% an imagesc figure that collapses across conditions to show LLE vs delay
% with MI magnitude.

close all;
clear all;
clc;

%% --- Script Configuration ---
counts_or_pdf = 'pdf'; % 'counts' or 'pdf'
invisible_y_axis = true;

% LLE histogram parameters
lle_range = [-1.5, 1.5];
n_bins_lle = 25;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins_lle), inf];

% LLE imagesc parameters
imagesc_lle_range = [-15 2]; % LLE range for the imagesc plot
n_bins_imagesc_lle_range = diff(imagesc_lle_range)+0;

% Mean rate histogram parameters
rate_range = [0, 5];
n_bins_rate = 25;
rate_bins = [linspace(rate_range(1), rate_range(2), n_bins_rate + 1), inf];

% MI parameters
mi_delay_for_histogram_samples = 200; % delay in samples for histogram column
delays_for_swarm_samples = [500]; % delays in samples for swarm plot
mi_range = [0, 3.5];
n_bins_mi = 25;
mi_bins = [linspace(mi_range(1), mi_range(2), n_bins_mi + 1), inf];

% LLE vs MI slice plot parameters
mi_delay_window_for_slice_samples = [25 100]; % delay window in samples

% Figure bin styling
outer_bin_width_multiplier = 1;
middle_ticks_lle = [0];
middle_ticks_rate = [0];
middle_ticks_mi = [0];

% Readable titles for conditions
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

%% Select experiment folder
fprintf('Please select the main experiment directory to plot...\n');
base_dir = uigetdir(pwd, 'Select Experiment Directory');
if isequal(base_dir, 0)
    fprintf('User cancelled. Exiting.\n');
    return;
end
fprintf('Selected directory: %s\n', base_dir);

%% Load data files
conditions_to_plot = {'no_adaptation', 'std_only', 'sfa_only', 'sfa_and_std'};
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
        fprintf('Error: Data file not found for condition "%s". Expected at: %s\n', condition_name, data_file);
        found_all = false;
    end
end
if ~found_all
    fprintf('Could not find all required data files. Exiting.\n');
    return;
end

%% Extract values
extracted_values = struct();
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    results = all_data.(condition_name).results;

    lles = [];
    mean_rates = [];
    mi_for_hist = [];
    mi_for_swarm = [];
    all_mi_vectors = {};
    delay_vec_stored = [];
    fs_stored = [];

    num_results = length(results);
    for j = 1:num_results
        res = results{j};
        if isstruct(res) && isfield(res, 'success') && res.success
            if isempty(fs_stored) && isfield(res, 'parameters') && isfield(res.parameters, 'fs')
                fs_stored = res.parameters.fs;
            end
            if isfield(res.simulation_output, 'LLE')
                lles(end+1) = res.simulation_output.LLE; %#ok<AGROW>
            end
            if isfield(res.simulation_output, 'mean_rate')
                mean_rates(end+1) = res.simulation_output.mean_rate; %#ok<AGROW>
            end
            % MI fields
            has_mi_fields = isfield(res.simulation_output, 'MI_vs_delay') && isfield(res.simulation_output, 'delay_vec');
            is_valid_mi = false;
            if has_mi_fields
                delay_vec = res.simulation_output.delay_vec;
                mi_vec = res.simulation_output.MI_vs_delay;
                if isvector(mi_vec) && ~any(isnan(mi_vec)) && isvector(delay_vec) && ~any(isnan(delay_vec)) && (length(mi_vec) == length(delay_vec))
                    is_valid_mi = true;
                end
            end
            if is_valid_mi
                all_mi_vectors{end+1} = mi_vec; %#ok<AGROW>
                % Extract MI at configured delay
                [~, delay_idx] = min(abs(delay_vec - mi_delay_for_histogram_samples));
                if ~isempty(delay_idx)
                    mi_for_hist(end+1) = mi_vec(delay_idx); %#ok<AGROW>
                end

                % Extract MI for swarm plot
                mi_at_delays = zeros(1, length(delays_for_swarm_samples));
                for k_delay = 1:length(delays_for_swarm_samples)
                    [~, idx] = min(abs(delay_vec - delays_for_swarm_samples(k_delay)));
                    mi_at_delays(k_delay) = mi_vec(idx);
                end
                mi_for_swarm(end+1, :) = mi_at_delays; %#ok<AGROW>

                if isempty(delay_vec_stored)
                    delay_vec_stored = delay_vec;
                end
            else
                if ~isempty(delay_vec_stored)
                    all_mi_vectors{end+1} = zeros(size(delay_vec_stored)); %#ok<AGROW>
                    mi_for_hist(end+1) = 0; %#ok<AGROW>
                    mi_for_swarm(end+1, :) = zeros(1, length(delays_for_swarm_samples)); %#ok<AGROW>
                end
            end
        end
    end
    extracted_values.(condition_name).lles = lles;
    extracted_values.(condition_name).mean_rates = mean_rates;
    extracted_values.(condition_name).mi_for_hist = mi_for_hist;
    extracted_values.(condition_name).mi_for_swarm = mi_for_swarm;
    extracted_values.(condition_name).all_mi_vectors = all_mi_vectors;
    extracted_values.(condition_name).delay_vec = delay_vec_stored;
    extracted_values.(condition_name).fs = fs_stored;
    fprintf('Condition "%s": %d LLEs, %d mean rates, %d MIs for histogram, %d full MI vectors.\n', ...
        condition_name, length(lles), length(mean_rates), length(mi_for_hist), length(all_mi_vectors));
end

%% Create main distributions figure (LLE and rate)
fig_main = figure('Name', 'LLE and Rate Distributions (Paired-Pulse)', 'Position', [100, 100, 650, 1020]);

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
    text(-0.2, 0.5, sprintf('%s', condition_titles(condition_name)), 'Units','normalized', 'Rotation',90, 'FontSize',22, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
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
% Build a pooled matrix with rows binned by LLE and columns by delay index; values aggregate median MI
pooled_delay = [];
pooled_mi = [];
pooled_lle = [];
template_delay_vec = [];
template_fs = [];

for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    delay_vec = extracted_values.(condition_name).delay_vec;
    fs = extracted_values.(condition_name).fs;
    mi_vectors = extracted_values.(condition_name).all_mi_vectors;
    lle_list = extracted_values.(condition_name).lles;
    % Align template
    if isempty(template_delay_vec) && ~isempty(delay_vec)
        template_delay_vec = delay_vec;
        template_fs = fs;
    end
    if isempty(delay_vec) || isempty(mi_vectors)
        continue
    end
    % Accumulate per-run MI vectors and corresponding LLEs (use nearest by index)
    for k = 1:length(mi_vectors)
        mi_vec = mi_vectors{k};
        if length(mi_vec) ~= length(delay_vec)
            continue
        end
        pooled_mi(end+1, :) = mi_vec(:)'; %#ok<AGROW>
        % Best effort: associate each MI vector with the corresponding LLE if available
        if ~isempty(lle_list)
            lle_val = lle_list(min(k, length(lle_list)));
        else
            lle_val = NaN;
        end
        pooled_lle(end+1, 1) = lle_val; %#ok<AGROW>
    end
end

fig_img = figure('Name', 'MI vs Delay by LLE (Collapsed Across Conditions)', 'Position', [643   321   730   580]);
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
    title('Mean paired-pulse MI vs delay and LLE')
    colorbar; ylabel(colorbar, 'MI (bits)');
end

%% LLE vs MI slice figure
fig_lle_mi_slice = figure('Name', 'LLE vs MI Slice', 'Position', [700, 400, 560, 420]);
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
        mi_valid = mi_averaged(valid_indices);

        % Use the same LLE bins as the imagesc plot for consistency
        lle_edges_slice = linspace(imagesc_lle_range(1), imagesc_lle_range(2), n_bins_imagesc_lle_range);
        n_lle_bins_slice = length(lle_edges_slice) - 1;
        
        lle_bin_centers = (lle_edges_slice(1:end-1) + lle_edges_slice(2:end)) / 2;
        mi_mean_in_bin = NaN(1, n_lle_bins_slice);
        mi_std_in_bin = NaN(1, n_lle_bins_slice);
        
        for b = 1:n_lle_bins_slice
            in_bin_indices = find(lle_valid >= lle_edges_slice(b) & lle_valid < lle_edges_slice(b+1));
            if ~isempty(in_bin_indices)
                mi_in_bin = mi_valid(in_bin_indices);
                mi_mean_in_bin(b) = mean(mi_in_bin, 'omitnan');
                mi_std_in_bin(b) = std(mi_in_bin, 'omitnan');
            end
        end
        
        % Plotting
        errorbar(lle_bin_centers, mi_mean_in_bin, mi_std_in_bin, 'o-', 'LineWidth', 1.5, 'CapSize', 4);
        xlabel('LLE (\lambda_1)');
        ylabel('Mutual Information (bits)');
        % title(sprintf('MI averaged over delays [%d, %d] samples', mi_delay_window_for_slice_samples(1), mi_delay_window_for_slice_samples(2)));
        box off;
        set(gca, 'FontSize', 14);
    end
end

%% Create Swarm Plot figure for MI at different delays
fig_swarm = figure('Name', 'MI Swarm Plot by Condition and Delay', 'Position', [200, 200, 1200, 600]);
tiledlayout(1, num_conditions, 'TileSpacing', 'compact');
all_ax_swarm = [];
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    mi_data = extracted_values.(condition_name).mi_for_swarm;

    if isempty(mi_data)
        nexttile;
        axis off;
        title(condition_titles(condition_name));
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        continue;
    end

    ax_swarm = nexttile;
    all_ax_swarm = [all_ax_swarm, ax_swarm]; %#ok<AGROW>

    % Prepare data for swarmchart
    num_points = size(mi_data, 1);
    num_delays = length(delays_for_swarm_samples);
    
    mi_flat = mi_data(:);
    delay_labels_categorical = categorical(repmat(delays_for_swarm_samples', num_points, 1));
    
    swarmchart(ax_swarm, delay_labels_categorical, mi_flat, 'filled');
    
    title(condition_titles(condition_name));
    xlabel('Delay (samples)');
    if i == 1
        ylabel('Mutual Information (bits)');
    else
        yticklabels({}); % Hide y-axis labels for other subplots
    end

end
if ~isempty(all_ax_swarm)
    linkaxes(all_ax_swarm, 'y');
end


%% Create Violin Plot figure for MI at different delays (using violinPlots2)
fig_violin = figure('Name', 'MI Violin Plot by Condition and Delay', 'Position', [200, 200, 1200, 600]);
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

    ax_violin = nexttile;
    all_ax_violin = [all_ax_violin, ax_violin]; %#ok<AGROW>

    if isempty(mi_data)
        axis off;
        title(condition_titles(condition_name));
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        continue;
    end
    
    % Use the custom violinplot, turn off its default median marker and use a shadow for the IQR
    delay_labels = string(delays_for_swarm_samples);
    violins = violinplot(mi_data, delay_labels, 'Parent', ax_violin, ...
        'ShowMedian', false, ...
        'QuartileStyle', 'shadow', ...
        'ShowBox', false, ...
        'ShowWhiskers', false);
    
    hold(ax_violin, 'on');

    % Calculate medians for each delay category
    medians = median(mi_data, 1, 'omitnan');
    
    % Loop through each violin object to draw a median line spanning its width
    for k = 1:length(violins)
        v = violins(k);
        y_median = medians(k);
        
        if isempty(v.ViolinPlot) || ~isvalid(v.ViolinPlot)
            continue;
        end
        
        violin_patch = v.ViolinPlot;
        x_data = violin_patch.XData;
        y_data = violin_patch.YData;
        
        % The YData of a patch goes from min to max, then max to min.
        % We can find the width by interpolating on each side.
        [y_max, max_idx] = max(y_data);
        
        if y_median >= min(y_data) && y_median <= y_max
            % Split data into the 'up' and 'down' side of the violin's contour
            y_up = y_data(1:max_idx);
            x_up = x_data(1:max_idx);
            y_down = y_data(max_idx:end);
            x_down = x_data(max_idx:end);

            % To use interp1, y vectors must be monotonic. Remove duplicates.
            [y_up_unique, ia_up] = unique(y_up);
            x_up_unique = x_up(ia_up);

            [y_down_unique, ia_down] = unique(y_down);
            x_down_unique = x_down(ia_down);
            
            % Interpolate to find the x-coordinates at the median y-value
            x1 = interp1(y_up_unique, x_up_unique, y_median, 'linear');
            x2 = interp1(y_down_unique, x_down_unique, y_median, 'linear');
    
            if isfinite(x1) && isfinite(x2)
                plot(ax_violin, [x1, x2], [y_median, y_median], 'k-', 'LineWidth', 2);
            end
        end
        
        % Manually add jitter to points at y=0, as the library code collapses them.
        if v.ShowData
            zero_data_indices = (mi_data(:, k) == 0);
            num_zero_points = sum(zero_data_indices);
        
            if num_zero_points > 0 && (0 >= min(y_data) && 0 <= y_max)
                % Re-use the unique-d x/y data from the median calculation above
                x1_at_zero = interp1(y_up_unique, x_up_unique, 0, 'linear', NaN);
                x2_at_zero = interp1(y_down_unique, x_down_unique, 0, 'linear', NaN);
        
                if isfinite(x1_at_zero) && isfinite(x2_at_zero)
                    width_at_zero = abs(x2_at_zero - x1_at_zero);
                    
                    % Generate jittered x-positions centered at the violin's position (k)
                    % jitter = (rand(num_zero_points, 1) - 0.5) * width_at_zero;
                    jitter = 2*(rand(num_zero_points, 1) - 0.5) * width_at_zero;
                    x_zeros = k + jitter;
                    
                    % Plot the points using properties from the original scatter plot
                    % to ensure they match in style.
                    if ~isempty(v.ScatterPlot) && isvalid(v.ScatterPlot)
                        scatter(ax_violin, x_zeros, zeros(num_zero_points, 1), ...
                                v.ScatterPlot.SizeData, v.ScatterPlot.CData, 'filled', ...
                                'MarkerFaceAlpha', v.ScatterPlot.MarkerFaceAlpha);
                    end
                end
            end
        end
    end
    hold(ax_violin, 'off');
    
    box(ax_violin, 'off'); % Turn off the box for all subplots
    
    % Use condition names as xlabels instead of titles, and remove x-ticks
    xlabel(ax_violin, condition_titles(condition_name));
    set(ax_violin, 'XTick', []);
    
    if i == 1
        ylabel('Mutual Information (bits)');
        yticks(ax_violin, [0 1 2 3]);
    else
        % Make y-axis invisible for subplots 2, 3, and 4
        set(ax_violin, 'YTick', [], 'YColor', 'none');
        ylabel('');
    end
end
if ~isempty(all_ax_violin)
    linkaxes(all_ax_violin, 'y');
    ylim(all_ax_violin, [0 3.2]); % Set matching y-limits for all plots
end


%% Save figures
output_dir_for_figs = fullfile(base_dir, 'analysis_plots');
if ~exist(output_dir_for_figs, 'dir')
    mkdir(output_dir_for_figs);
end

% fprintf('\nSaving figures to: %s\n', output_dir_for_figs);
save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_LLE_rate_MI_distributions_v1', fig_main.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_MI_vs_delay_LLE_imagesc_v1', fig_img.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_LLE_vs_MI_slice_v1', fig_lle_mi_slice.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_MI_swarm_by_condition_v1', fig_swarm.Number, {'fig', 'svg', 'png'});
save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_MI_violin_by_condition_v1', fig_violin.Number, {'fig', 'svg', 'png'});


fprintf('Plotting complete.\n');


