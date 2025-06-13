function fig_handle = sensitivity_plot(param_file, hist_edges_with_inf)
    %SENSITIVITY_PLOT Creates visualization of LLE sensitivity analysis results
    %
    % Inputs:
    %   param_file           - path to the parameter-specific .mat file
    %   hist_edges_with_inf  - vector of bin edges for LLE histograms (can include -inf, inf)
    %
    % Outputs:
    %   fig_handle           - handle to the generated (invisible) figure
    
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
    all_lle_values_by_level = cell(n_levels, 1);
    
    for level_idx = 1:n_levels
        lle_values_level = [];
        
        for rep_idx = 1:n_reps
            result = results{level_idx, rep_idx};
            
            if isfield(result, 'success') && result.success && isfield(result, 'LLE')
                if isnan(result.LLE)
                    lle_values_level(end+1) = 1e3; % Assign a large value for NaN LLEs
                else
                    lle_values_level(end+1) = result.LLE;
                end
                success_counts(level_idx) = success_counts(level_idx) + 1;
            end
        end
        
        all_lle_values_by_level{level_idx} = lle_values_level;
        
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
    if ~isempty(lle_values_all)
        fprintf('LLE range: [%.3f, %.3f]\n', min(lle_values_all), max(lle_values_all));
    else
        fprintf('No successful runs with LLE values found.\n');
    end
    
    % Prepare finite y-coordinates for imagesc
    internal_finite_edges = hist_edges_with_inf(2:end-1);
    y_coords_for_plot = zeros(num_hist_bins, 1);

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
    fig_handle = figure('Position', [100, 100, 600, 600], 'Visible', 'off');
    imagesc(param_levels, y_coords_for_plot, lle_histogram_matrix);
    % colorbar;
    caxis([0 n_reps]); % Set colormap extents from 0 to n_reps
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel('Largest Lyapunov Exponent (LLE)');
    % title(sprintf('LLE Distribution vs %s', strrep(param_name, '_', '\\_')));
    axis xy; % Flip y-axis so smaller LLE values are at bottom
    % colormap(hot); % Use blue-yellow colormap
    colormap(parula)

    % Add grid lines for better readability
    % hold on;
    % for i = 1:length(param_levels)
    %     line([param_levels(i), param_levels(i)], [plot_y_min, plot_y_max], ...
    %          'Color', 'k', 'LineWidth', 0.2);
    % end
    % hold off;
    drawnow; % Force update for the plot

    % Save the figure
    % if nargin >= 3 && ~isempty(output_dir)
    %     fig_filename = fullfile(output_dir, sprintf('sensitivity_%s.png', param_name));
    %     saveas(gcf, fig_filename);
    %     fig_filename_fig = fullfile(output_dir, sprintf('sensitivity_%s.fig', param_name));
    %     saveas(gcf, fig_filename_fig);
    %     fprintf('Saved heatmap: %s\n', fig_filename);
    % end
end 