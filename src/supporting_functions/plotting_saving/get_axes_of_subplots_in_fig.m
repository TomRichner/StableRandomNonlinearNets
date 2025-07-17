function axes_in_fig = get_axes_of_subplots_in_fig(fig_handle, sort_direction)
% get_axes_of_subplots_in_fig - Finds and sorts all main subplot axes in a figure.
%
% This function finds all axes objects in a given figure, but excludes
% common non-subplot axes like legends, colorbars, and helper axes. It then
% sorts the axes by their position on the figure.
%
% Syntax:
%   axes_in_fig = get_axes_of_subplots_in_fig(fig_handle, sort_direction)
%
% Inputs:
%   fig_handle - (Optional) Handle to the figure object. Defaults to gcf.
%   sort_direction - (Optional) 'LeftRight' (default) or 'TopDown'.
%
% Outputs:
%   axes_in_fig - A sorted array of handles to the identified subplot axes.

if nargin < 1
    fig_handle = gcf;
end
if nargin < 2
    sort_direction = 'LeftRight'; % Default sort order
end

% Find all axes, excluding some by tag
lAxes_unsorted = findall(fig_handle, 'type', 'axes', ...
    '-not', {'Tag', 'legend'}, ...
    '-not', {'Tag', 'Colorbar'}, ...
    '-not', {'Tag', 'scalebar_text_ax'});

if isempty(lAxes_unsorted)
    axes_in_fig = lAxes_unsorted;
    return;
end

% Sort axes by position
Pos = NaN(length(lAxes_unsorted), 2);
for ia = 1:length(lAxes_unsorted)
    ax = lAxes_unsorted(ia);
    Pos(ia, :) = [ax.InnerPosition(1),  ax.InnerPosition(2)];
end

% Invert Y position for sorting because lower Y is higher on the figure
Pos(:, 2) = -Pos(:, 2); 

switch sort_direction
    case 'LeftRight'
        [~, ind] = sortrows(Pos, [2, 1]); % Sort by Y (top-to-bottom), then X (left-to-right)
    case 'TopDown'
        [~, ind] = sortrows(Pos, [1, 2]); % Sort by X (left-to-right), then Y (top-to-bottom)
    otherwise
        warning('Unknown sort direction: %s. Using default ''LeftRight''.', sort_direction);
        [~, ind] = sortrows(Pos, [2, 1]);
end

axes_in_fig = lAxes_unsorted(ind);

end
