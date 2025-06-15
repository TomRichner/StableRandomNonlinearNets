function fig_handle = sensitivity_plot(param_file, hist_edges_with_inf, variable_to_plot, y_axis_label)
    %SENSITIVITY_PLOT Creates visualization of sensitivity analysis results
    %
    % Inputs:
    %   param_file           - path to the parameter-specific .mat file
    %   hist_edges_with_inf  - vector of bin edges for histograms (can include -inf, inf)
    %   variable_to_plot     - name of the field to plot from the results (e.g., 'LLE', 'mean_rate')
    %   y_axis_label         - label for the y-axis
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
    
    % Extract values
    num_hist_bins = length(hist_edges_with_inf) - 1;
    histogram_matrix = zeros(num_hist_bins, n_levels);
    values_all = [];
    success_counts = zeros(n_levels, 1);
    all_values_by_level = cell(n_levels, 1);
    
    for level_idx = 1:n_levels
        values_level = [];
        
        for rep_idx = 1:n_reps
            result = results{level_idx, rep_idx};
            
            if isfield(result, 'success') && result.success && isfield(result, variable_to_plot)
                val = result.(variable_to_plot);
                if isnan(val)
                    values_level(end+1) = 1e3; % Assign a large value for NaN values
                else
                    values_level(end+1) = val;
                end
                success_counts(level_idx) = success_counts(level_idx) + 1;
            end
        end
        
        all_values_by_level{level_idx} = values_level;
        
        % Create histogram for this level
        if ~isempty(values_level)
            [counts, ~] = histcounts(values_level, hist_edges_with_inf);
            histogram_matrix(:, level_idx) = counts;
            values_all = [values_all, values_level];
        end
    end
    
    % Calculate success rate
    total_success = sum(success_counts);
    total_attempts = n_levels * n_reps;
    success_rate = total_success / total_attempts;
    
    fprintf('Overall success rate for %s: %d/%d (%.1f%%)\n', variable_to_plot, total_success, total_attempts, 100*success_rate);
    if ~isempty(values_all)
        fprintf('%s range: [%.3f, %.3f]\n', variable_to_plot, min(values_all), max(values_all));
    else
        fprintf('No successful runs with %s values found.\n', variable_to_plot);
    end

    % Extract finite range from hist_edges_with_inf for custom tick labels
    finite_edges = hist_edges_with_inf(~isinf(hist_edges_with_inf));
    if length(finite_edges) >= 2
        value_min = min(finite_edges);
        value_max = max(finite_edges);
    else
        % Fallback if no finite edges found
        value_min = -1;
        value_max = 1;
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
    imagesc(param_levels, y_coords_for_plot, histogram_matrix);
    % colorbar;
    caxis([0 n_reps]); % Set colormap extents from 0 to n_reps
    xlabel(sprintf('%s', strrep(param_name, '_', '\\_')));
    ylabel(y_axis_label, 'Interpreter', 'latex', 'FontSize', 22);
    % title(sprintf('Distribution of %s vs %s', variable_to_plot, strrep(param_name, '_', '\\_')));
    axis xy; % Flip y-axis so smaller values are at bottom
    % colormap(hot);
    colormap(parula)

    % Add custom y-tick labels with a fixed set of intermediate values, plus
    % labels for the outermost bins that collect outlier values.
    if exist('step_size', 'var')
        % Define desired intermediate ticks.
        intermediate_ticks = unique([-1, 0, 1, round(value_min), round(value_max)]);
        if strcmp(variable_to_plot, 'mean_rate')
            intermediate_ticks = unique([0, 10, 20, 30, 40]);
        end
        
        % Filter out any intermediate ticks that are outside the finite range.
        intermediate_ticks = intermediate_ticks(intermediate_ticks >= value_min & intermediate_ticks <= value_max);

        % If an intermediate tick is at the minimum value, remove it to avoid
        % a redundant tick, as the first bin's tick already represents this.
        if ~isempty(intermediate_ticks)
            intermediate_ticks(abs(intermediate_ticks - value_min) < 1e-9) = [];
        end

        % Define the final ticks: user-specified intermediate ticks plus ticks for
        % the top and bottom bins (which collect outliers).
        % The y-coordinates for the top/bottom bins are their calculated centers,
        % while intermediate ticks are placed at their own value on the y-axis.
        final_yticks = unique([y_coords_for_plot(1); ...
                               intermediate_ticks(:); ...
                               y_coords_for_plot(end)]);

        % Generate labels for the final ticks.
        final_ylabels = cell(size(final_yticks));
        for k = 1:length(final_yticks)
            tick_val = final_yticks(k);
            % Use a small tolerance for float comparison.
            if abs(tick_val - y_coords_for_plot(1)) < 1e-6
                if value_min == 0
                    final_ylabels{k} = sprintf('%.1f', value_min);
                else
                    final_ylabels{k} = sprintf('≤%.1f', value_min);
                end
            elseif abs(tick_val - y_coords_for_plot(end)) < 1e-6
                final_ylabels{k} = sprintf('≥%.1f', value_max);
            else
                % For intermediate ticks, the label is the value itself.
                final_ylabels{k} = sprintf('%.1f', tick_val);
            end
        end

        % Apply the new ticks and labels.
        yticks(final_yticks);
        yticklabels(final_ylabels);
    end

    % Ensure the first parameter level is always an x-tick
    current_xticks = xticks;
    if ~isempty(param_levels)
        if ~any(abs(current_xticks - param_levels(1)) < 1e-9) % Use a tolerance for float comparison
            xticks(sort([param_levels(1), current_xticks]));
        end
    end

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