function sensitivity_plot(param_file, hist_edges_with_inf, output_dir)
    %SENSITIVITY_PLOT Creates visualization of LLE sensitivity analysis results
    %
    % Inputs:
    %   param_file           - path to the parameter-specific .mat file
    %   hist_edges_with_inf  - vector of bin edges for LLE histograms (can include -inf, inf)
    %   output_dir           - directory to save plots
    
    % Load the parameter results
    fprintf('Loading: %s\n', param_file);
    data = load(param_file);
    
    if ~isfield(data, 'results_reshaped') || ~isfield(data, 'metadata')
        error('Expected fields "results_reshaped" and "metadata" not found in %s', param_file);
    end
    
    results = data.results_reshaped;
    metadata = data.metadata;
    
    param_name = metadata.param_name;
    param_levels = metadata.param_levels;
    n_levels = metadata.n_levels;
    n_reps = metadata.n_reps;
    
    fprintf('Processing %s: %d levels, %d reps per level\n', param_name, n_levels, n_reps);
    
    % Extract LLE values
    num_hist_bins = length(hist_edges_with_inf) - 1;
    lle_histogram_matrix = zeros(num_hist_bins, n_levels);
    lle_values_all = [];
    success_counts = zeros(n_levels, 1);
    
    for level_idx = 1:n_levels
        lle_values_level = [];
        
        for rep_idx = 1:n_reps
            result = results{level_idx, rep_idx};
            
            if isfield(result, 'success') && result.success && isfield(result, 'LLE')
                lle_values_level(end+1) = result.LLE;
                success_counts(level_idx) = success_counts(level_idx) + 1;
            end
        end
        
        % Create histogram for this level
        if ~isempty(lle_values_level)
            [counts, ~] = histcounts(lle_values_level, hist_edges_with_inf);
            lle_histogram_matrix(:, level_idx) = counts;
            lle_values_all = [lle_values_all, lle_values_level];
        end
    end
    
    % Calculate success rate
    total_success = sum(success_counts);
    total_attempts = n_levels * n_reps;
    success_rate = total_success / total_attempts;
    
    fprintf('Overall success rate: %d/%d (%.1f%%)\n', total_success, total_attempts, 100*success_rate);
    fprintf('LLE range: [%.3f, %.3f]\n', min(lle_values_all), max(lle_values_all));
    
    % Prepare finite y-coordinates for imagesc
    internal_finite_edges = hist_edges_with_inf(2:end-1);
    y_coords_for_plot = zeros(num_hist_bins, 1);
    plot_y_min = 0; % Default, will be updated
    plot_y_max = 1; % Default, will be updated

    if num_hist_bins == 0
        % Should not happen with current caller settings
        warning('No histogram bins found.');
        % Assign some defaults to prevent further errors, though plot will be meaningless
        y_coords_for_plot = [0];
        plot_y_min = -1;
        plot_y_max = 1;
    elseif isempty(internal_finite_edges) % Only one bin: (-inf, inf)
        % This implies hist_edges_with_inf was [-inf, inf]
        % num_hist_bins will be 1
        y_coords_for_plot(1) = 0; % Arbitrary center for the single bin
        plot_y_min = -1; % Arbitrary range for plotting
        plot_y_max = 1;
    elseif length(internal_finite_edges) == 1 % Two bins: (-inf, edge] and (edge, inf)
        % num_hist_bins will be 2
        step_size = abs(internal_finite_edges(1) * 0.2);
        if step_size == 0, step_size = 1.0; end % Avoid zero step_size
        
        y_coords_for_plot(1) = internal_finite_edges(1) - step_size/2;
        y_coords_for_plot(2) = internal_finite_edges(1) + step_size/2;
        plot_y_min = internal_finite_edges(1) - step_size;
        plot_y_max = internal_finite_edges(1) + step_size;
    else % Standard case: >=3 bins, including at least one internal bin
        step_size = internal_finite_edges(2) - internal_finite_edges(1);
        if step_size <= 0 % Should not happen with linspace if range is valid
            warning('Calculated step_size is not positive. Using default.');
            step_size = abs(mean(diff(internal_finite_edges))); % Try average diff
            if step_size <= 0 || isnan(step_size) % Fallback
                 step_size = 1.0;
            end
        end

        y_coords_for_plot(1) = internal_finite_edges(1) - step_size/2;
        for k = 1:length(internal_finite_edges)-1
            y_coords_for_plot(k+1) = (internal_finite_edges(k) + internal_finite_edges(k+1))/2;
        end
        y_coords_for_plot(end) = internal_finite_edges(end) + step_size/2;
        
        plot_y_min = internal_finite_edges(1) - step_size;
        plot_y_max = internal_finite_edges(end) + step_size;
    end

    % Create the main visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Main imagesc plot (occupies a 2x2 block in a 2x3 grid: left 2/3 of figure)
    subplot(2, 3, [1, 2, 4, 5]); % Spans cells 1,2 (top row, cols 1-2) and 4,5 (bottom row, cols 1-2)
    imagesc(param_levels, y_coords_for_plot, lle_histogram_matrix);
    colorbar;
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel('Largest Lyapunov Exponent (LLE)');
    title(sprintf('LLE Distribution vs %s', strrep(param_name, '_', '\\_')));
    axis xy; % Flip y-axis so smaller LLE values are at bottom
    colormap(hot);
    
    % Add grid lines for better readability
    hold on;
    for i = 1:length(param_levels)
        line([param_levels(i), param_levels(i)], [plot_y_min, plot_y_max], ...
             'Color', 'k', 'LineWidth', 0.2, 'Alpha', 0.3);
    end
    hold off;
    drawnow; % Force update for the first subplot

    % Debugging output for success rate plot
    fprintf('Debug MATPLOT: Param: %s. Attempting to plot success rate.\\n', param_name);
    fprintf('Debug MATPLOT: Num param_levels: %d. Num success_counts: %d. n_reps: %d\\n', length(param_levels), length(success_counts), n_reps);
    if ~isempty(param_levels) && length(param_levels) > 0
        fprintf('Debug MATPLOT: param_levels(1) = %f\\n', param_levels(1));
    else
        fprintf('Debug MATPLOT: param_levels is empty or not populated.\\n');
    end
    if ~isempty(success_counts) && length(success_counts) > 0
        fprintf('Debug MATPLOT: success_counts(1) = %d\\n', success_counts(1));
    else
        fprintf('Debug MATPLOT: success_counts is empty or not populated.\\n');
    end
    
    % Success rate plot (DIAGNOSTIC - simplified)
    subplot(2, 3, 3);
    if ~isempty(param_levels) && ~isempty(success_counts) && length(param_levels) == length(success_counts)
        plot(param_levels, success_counts / n_reps * 100, 'o-');
        ylabel('Success Rate (%)');
        title(sprintf('Success Rate (Test Plot) for %s', strrep(param_name, '_', '\\_')));
    else
        text(0.5,0.5, 'Data missing for success rate', 'HorizontalAlignment', 'center');
        title(sprintf('Success Rate (Data Missing) for %s', strrep(param_name, '_', '\\_')));
    end
    ylim([0, 105]); % Ensure 100% is visible if plotted
    grid on;
    drawnow; % Force update

    % Debugging output for LLE distribution plot
    fprintf('Debug MATPLOT: Param: %s. Attempting to plot overall LLE distribution.\\n', param_name);
    fprintf('Debug MATPLOT: Num lle_values_all: %d. Num lle_bins: %d\\n', length(lle_values_all), length(hist_edges_with_inf));
    if ~isempty(lle_values_all) && length(lle_values_all) > 0
        fprintf('Debug MATPLOT: lle_values_all(1) = %f\\n', lle_values_all(1));
    else
        fprintf('Debug MATPLOT: lle_values_all is empty or not populated for histogram.\\n');
    end

    % Overall LLE distribution (DIAGNOSTIC - simplified)
    subplot(2, 3, 6);
    if ~isempty(lle_values_all)
        % Using a simple plot instead of histogram for diagnostics
        plot(sort(lle_values_all), 'o-'); 
        title(sprintf('Sorted LLE Values (Test Plot) (Count: %d)', length(lle_values_all)));
        ylabel('LLE Value');
        xlabel('Sorted Index');
    else
        text(0.5,0.5, 'lle_values_all is empty', 'HorizontalAlignment', 'center');
        title('Overall LLE Distribution (Data Missing)');
    end
    grid on;
    % Add vertical line at LLE = 0 if data exists
    if ~isempty(lle_values_all)
        hold on;
        % Find index near LLE=0 if plotting sorted LLEs for context, or just plot yline
        yline(0, 'r--', 'LineWidth', 2);
        hold off;
    end
    drawnow; % Force update
    
    % Add text annotations
    chaos_fraction = sum(lle_values_all > 0) / length(lle_values_all);
    text(0.05, 0.95, sprintf('Chaos fraction: %.2f', chaos_fraction), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    
    % Overall title
    condition_str = '';
    if isfield(metadata, 'condition') && isfield(metadata.condition, 'name')
        condition_str_raw = strrep(metadata.condition.name, '_', ' ');
        condition_str = sprintf(' (Condition: %s)', condition_str_raw);
    end
    
    sgtitle(sprintf('Sensitivity Analysis: %s%s (n=%d successful runs)', ...
                   strrep(param_name, '_', '\\_'), condition_str, total_success), 'FontSize', 14);
    
    % Save the figure
    if nargin >= 3 && ~isempty(output_dir)
        fig_filename = fullfile(output_dir, sprintf('sensitivity_%s.png', param_name));
        saveas(gcf, fig_filename);
        
        fig_filename_fig = fullfile(output_dir, sprintf('sensitivity_%s.fig', param_name));
        saveas(gcf, fig_filename_fig);
        
        fprintf('Saved plots: %s\n', fig_filename);
    end
    
    % Create a detailed statistics plot
    figure('Position', [150, 150, 1000, 600]);
    
    % Mean LLE vs parameter level
    subplot(2, 2, 1);
    mean_lle = zeros(n_levels, 1);
    std_lle = zeros(n_levels, 1);
    
    for level_idx = 1:n_levels
        lle_values_level = [];
        for rep_idx = 1:n_reps
            result = results{level_idx, rep_idx};
            if isfield(result, 'success') && result.success && isfield(result, 'LLE')
                lle_values_level(end+1) = result.LLE;
            end
        end
        
        if ~isempty(lle_values_level)
            mean_lle(level_idx) = mean(lle_values_level);
            std_lle(level_idx) = std(lle_values_level);
        else
            mean_lle(level_idx) = NaN;
            std_lle(level_idx) = NaN;
        end
    end
    
    errorbar(param_levels, mean_lle, std_lle, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel('Mean LLE ± STD');
    title('Mean LLE vs Parameter');
    grid on;
    
    % Add horizontal line at LLE = 0
    hold on;
    line(xlim, [0, 0], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    hold off;
    
    % Chaos probability vs parameter level
    subplot(2, 2, 2);
    chaos_prob = zeros(n_levels, 1);
    
    for level_idx = 1:n_levels
        lle_values_level = [];
        for rep_idx = 1:n_reps
            result = results{level_idx, rep_idx};
            if isfield(result, 'success') && result.success && isfield(result, 'LLE')
                lle_values_level(end+1) = result.LLE;
            end
        end
        
        if ~isempty(lle_values_level)
            chaos_prob(level_idx) = sum(lle_values_level > 0) / length(lle_values_level);
        else
            chaos_prob(level_idx) = NaN;
        end
    end
    
    plot(param_levels, chaos_prob * 100, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel('Chaos Probability (%)');
    title('Probability of Chaotic Dynamics');
    ylim([0, 100]);
    grid on;
    
    % Parameter value distribution
    subplot(2, 2, 3);
    histogram(param_levels, min(20, n_levels), 'Normalization', 'count');
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel('Count');
    title('Parameter Sampling Distribution');
    grid on;
    
    % Summary statistics
    subplot(2, 2, 4);
    axis off;
    
    stats_text = {
        sprintf('Parameter: %s', param_name);
        sprintf('Range: [%.3f, %.3f]', min(param_levels), max(param_levels));
        sprintf('Levels: %d', n_levels);
        sprintf('Reps per level: %d', n_reps);
        sprintf('Total runs: %d', total_attempts);
        sprintf('Successful runs: %d (%.1f%%)', total_success, 100*success_rate);
        '';
        sprintf('LLE Statistics:');
        sprintf('  Mean: %.4f', nanmean(lle_values_all));
        sprintf('  Std: %.4f', nanstd(lle_values_all));
        sprintf('  Min: %.4f', min(lle_values_all));
        sprintf('  Max: %.4f', max(lle_values_all));
        sprintf('  Chaos fraction: %.3f', chaos_fraction);
    };
    
    text(0.1, 0.9, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', 10, 'FontFamily', 'monospace');
    
    sgtitle(sprintf('Detailed Statistics: %s', strrep(param_name, '_', '\\_')), 'FontSize', 14);
    
    % Save the detailed statistics figure
    if nargin >= 3 && ~isempty(output_dir)
        fig_filename_stats = fullfile(output_dir, sprintf('sensitivity_%s_stats.png', param_name));
        saveas(gcf, fig_filename_stats);
        
        fig_filename_stats_fig = fullfile(output_dir, sprintf('sensitivity_%s_stats.fig', param_name));
        saveas(gcf, fig_filename_stats_fig);
        
        fprintf('Saved detailed stats: %s\n', fig_filename_stats);
    end
    
    fprintf('Plotting complete for %s\n', param_name);
end 