% Figure defaults
set(groot, 'DefaultFigureRenderer', 'painters'); % painters is better for SVG
set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesLineWidth', 2);
set_lines_no_red_cmap % remove red from the default lines colormap because red is reserved for inhibitory neurons