% Add paths to necessary functions
addpath('../src/algorithms/Kaplan_Yorke');
addpath('../src/algorithms/Lyapunov');
addpath('../src/algorithms/analytic_fixed_point_spectrum');
addpath('../src/algorithms/ode_solvers');
addpath('../src/SRNN');
addpath('../src/SRNN/SRNN_utils');
addpath('../src/supporting_functions/plotting_saving');

% Figure defaults
set(groot, 'DefaultFigureRenderer', 'painters'); % painters is better for SVG
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesLineWidth', 2);
set_lines_no_red_cmap % remove red from the default lines colormap because red is reserved for inhibitory neurons
