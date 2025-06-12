addpath('supporting_functions');
addpath('ode_solvers');
addpath('fig_files')
addpath(['fig_files' filesep 'IEDs'])
addpath(['external_code' filesep 'interp1qr'])
addpath(['analytic_fixed_point_spectrum'])
set(groot, 'DefaultFigureRenderer', 'painters');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesLineWidth', 2);
set_lines_no_red_cmap
