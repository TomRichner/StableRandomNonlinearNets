addpath(fileparts(mfilename('fullpath'))); % add the directory that set_path_and_fig_defaults.m is in
addpath('supporting_functions');
addpath('ode_solvers');
addpath(['analytic_fixed_point_spectrum'])

set(groot, 'DefaultFigureRenderer', 'painters'); % painters is better for SVG
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesLineWidth', 2);
set_lines_no_red_cmap % remove red from the default lines colormap because red is reserved for inhibitory neurons
