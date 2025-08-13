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

% Mean rate histogram parameters
rate_range = [0, 10];
n_bins_rate = 25;
rate_bins = [linspace(rate_range(1), rate_range(2), n_bins_rate + 1), inf];

% MI parameters
mi_delay_for_histogram_samples = 200; % delay in samples for histogram column
mi_range = [0, 2.5];
n_bins_mi = 25;
mi_bins = [linspace(mi_range(1), mi_range(2), n_bins_mi + 1), inf];

% Figure bin styling
outer_bin_width_multiplier = 1;
middle_ticks_lle = [0];
middle_ticks_rate = [0, 5];
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
                if isempty(delay_vec_stored)
                    delay_vec_stored = delay_vec;
                end
            else
                if ~isempty(delay_vec_stored)
                    all_mi_vectors{end+1} = zeros(size(delay_vec_stored)); %#ok<AGROW>
                    mi_for_hist(end+1) = 0; %#ok<AGROW>
                end
            end
        end
    end
    extracted_values.(condition_name).lles = lles;
    extracted_values.(condition_name).mean_rates = mean_rates;
    extracted_values.(condition_name).mi_for_hist = mi_for_hist;
    extracted_values.(condition_name).all_mi_vectors = all_mi_vectors;
    extracted_values.(condition_name).delay_vec = delay_vec_stored;
    extracted_values.(condition_name).fs = fs_stored;
    fprintf('Condition "%s": %d LLEs, %d mean rates, %d MIs for histogram, %d full MI vectors.\n', ...
        condition_name, length(lles), length(mean_rates), length(mi_for_hist), length(all_mi_vectors));
end

%% Create main distributions figure (LLE, rate, MI)
fig_main = figure('Name', 'LLE, Rate, and MI Distributions (Paired-Pulse)', 'Position', [100, 100, 900, 1020]);

% Column 1: LLE
ax_lle = gobjects(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_lle(i) = subplot(num_conditions, 3, 3*i-2);
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
    ax_rate(i) = subplot(num_conditions, 3, 3*i-1);
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

% Column 3: MI at chosen delay
ax_mi = gobjects(num_conditions, 1);
for i = 1:num_conditions
    condition_name = conditions_to_plot{i};
    ax_mi(i) = subplot(num_conditions, 3, 3*i);
    [counts, ~] = histcounts(extracted_values.(condition_name).mi_for_hist, mi_bins);
    % Finite bins for MI (lower bound at zero)
    std_width_mi = (mi_range(2) - mi_range(1)) / n_bins_mi;
    middle_edges_mi = linspace(mi_range(1), mi_range(2), n_bins_mi + 1);
    first_edge = middle_edges_mi(1);
    last_edge = middle_edges_mi(end) + outer_bin_width_multiplier * std_width_mi;
    mi_bins_finite = [first_edge, middle_edges_mi(2:end), last_edge];
    if strcmpi(counts_or_pdf, 'pdf')
        plot_values = counts ./ (sum(counts) + eps);
    else
        plot_values = counts;
    end
    histogram('BinEdges', mi_bins_finite, 'BinCounts', plot_values, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    if i == num_conditions
        xlabel(sprintf('MI (bits)'),'FontSize',22);
    end
    if invisible_y_axis, set(gca, 'YTick', [], 'YTickLabel', [], 'YColor', 'none'); end
    box off;
    ax = gca; xlim([mi_bins_finite(1), mi_bins_finite(end)]);
    center_last_bin = (mi_bins_finite(end-1) + mi_bins_finite(end)) / 2;
    final_ticks = [middle_ticks_mi, center_last_bin];
    new_labels = [arrayfun(@(x) sprintf('%.1f', x), middle_ticks_mi, 'UniformOutput', false), {sprintf('>= %.1f', mi_range(2))}];
    xticks(ax, final_ticks); xticklabels(ax, new_labels); ax.XTickLabelRotation = 45;
end

% Link y axes
linkaxes([ax_lle; ax_rate; ax_mi], 'y');

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

fig_img = figure('Name', 'MI vs Delay by LLE (Collapsed Across Conditions)', 'Position', [150, 150, 900, 500]);
if isempty(pooled_mi) || isempty(template_delay_vec)
    axes; axis off; text(0.5,0.5,'No MI data available for imagesc','HorizontalAlignment','center');
else
    % Bin LLE into edges and compute median MI in each LLE bin per delay
    lle_edges = linspace(lle_range(1), lle_range(2), 15);
    n_lle_bins = length(lle_edges)-1;
    n_delays = length(template_delay_vec);
    mi_img = NaN(n_lle_bins, n_delays);

    for b = 1:n_lle_bins
        in_bin = pooled_lle >= lle_edges(b) & pooled_lle < lle_edges(b+1);
        if any(in_bin)
            mi_in_bin = pooled_mi(in_bin, :);
            mi_img(b, :) = median(mi_in_bin, 1, 'omitnan');
        end
    end

    imagesc((template_delay_vec./template_fs), (lle_edges(1:end-1)+lle_edges(2:end))/2, mi_img);
    axis xy; colormap(hot);
    xlabel('Delay (s)','FontSize',22);
    ylabel('LLE (binned center)','FontSize',22);
    title('Median MI vs Delay across LLE bins','FontSize',18);
    colorbar; ylabel(colorbar, 'MI (bits)');
end

%% Save figures
output_dir_for_figs = fullfile(base_dir, 'analysis_plots');
if ~exist(output_dir_for_figs, 'dir')
    mkdir(output_dir_for_figs);
end

fprintf('\nSaving figures to: %s\n', output_dir_for_figs);
if exist('save_some_figs_to_folder_2','file')
    save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_LLE_rate_MI_distributions_v1', fig_main.Number, {'fig', 'svg', 'png'});
    save_some_figs_to_folder_2(output_dir_for_figs, 'PPMI_MI_vs_delay_LLE_imagesc_v1', fig_img.Number, {'fig', 'svg', 'png'});
else
    % Fallback
    savefig(fig_main, fullfile(output_dir_for_figs, 'PPMI_LLE_rate_MI_distributions_v1.fig'));
    saveas(fig_main, fullfile(output_dir_for_figs, 'PPMI_LLE_rate_MI_distributions_v1.png'));
    savefig(fig_img, fullfile(output_dir_for_figs, 'PPMI_MI_vs_delay_LLE_imagesc_v1.fig'));
    saveas(fig_img, fullfile(output_dir_for_figs, 'PPMI_MI_vs_delay_LLE_imagesc_v1.png'));
end

fprintf('Plotting complete.\n');
beep; pause(0.5); beep;


