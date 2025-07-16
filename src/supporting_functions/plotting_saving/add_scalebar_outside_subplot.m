function add_scalebar_outside_subplot(ax, scale_bar_height_data, label_text)
% add_scalebar_outside_subplot - Adds a vertical scale bar with a label to the right of a subplot.
%
% Inputs:
%   ax - Handle to the axes object.
%   scale_bar_height_data - The height of the scale bar in data units.
%   label_text - The string to use as the label for the scale bar.

% --- Parameters for appearance ---
x_offset = 0.015;       % Horizontal offset from the right edge of the axes
lineWidth = 2.5;
y_start_data = 0;      % The y-position (in data units) where the scale bar starts
text_x_offset_from_bar = 0.00; % Horizontal offset of text from scale bar

% To place annotation outside axis, we need to convert from data to normalized figure units.
ax_pos = get(ax, 'Position'); % [left, bottom, width, height] in normalized figure units
ax_ylim = get(ax, 'YLim');
y_range = diff(ax_ylim);

if y_range == 0 % Avoid division by zero if ylim is flat
    return;
end

% We will place the bar from y_start_data to y_start_data + scale_bar_height_data, then convert to figure coordinates.
y_end_data = y_start_data + scale_bar_height_data;

% Convert y data coordinates to normalized figure coordinates for the annotation
y_start_fig = ax_pos(2) + ((y_start_data - ax_ylim(1)) / y_range) * ax_pos(4);
y_end_fig = ax_pos(2) + ((y_end_data - ax_ylim(1)) / y_range) * ax_pos(4);

% Position the scale bar to the right of the plot area
x_pos_fig = ax_pos(1) + ax_pos(3) + x_offset;

% Check if scale bar would be reasonably on screen.
% This is a heuristic check to avoid placing annotations way off figure.
% if y_start_fig > 0.05 && y_end_fig < 0.95
    % Draw the scale bar using annotation
    annotation('line', [x_pos_fig, x_pos_fig], [y_start_fig, y_end_fig], 'LineWidth', lineWidth, 'Color', 'k');

    % Add the rotated text label. Using text() on a temporary full-figure axis is a reliable way to get rotated text.
    text_x_fig = x_pos_fig + text_x_offset_from_bar;
    text_y_fig = (y_start_fig + y_end_fig) / 2;
    
    fig = get(ax, 'Parent');
    h_ax_for_text = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'Tag', 'scalebar_text_ax');
    text(h_ax_for_text, text_x_fig, text_y_fig, label_text, ...
        'Rotation', -90, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
    
    % Set current axes back to the main plot
    axes(ax);
% end

end
