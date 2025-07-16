function axes_in_fig = get_axes_of_subplots_in_fig(fig_handle)
% get_axes_of_subplots_in_fig - Finds all main subplot axes in a figure.
%
% This function finds all axes objects in a given figure, but excludes
% common non-subplot axes like legends, colorbars, and helper axes created
% by functions like `add_scalebar_outside_subplot`.
%
% Syntax:
%   axes_in_fig = get_axes_of_subplots_in_fig(fig_handle)
%
% Inputs:
%   fig_handle - (Optional) Handle to the figure object. Defaults to gcf.
%
% Outputs:
%   axes_in_fig - An array of handles to the identified subplot axes.

if nargin < 1
    fig_handle = gcf;
end

axes_in_fig = findall(fig_handle, 'type', 'axes', ...
    '-not', {'Tag', 'legend'}, ...
    '-not', {'Tag', 'Colorbar'}, ...
    '-not', {'Tag', 'scalebar_text_ax'});

end
