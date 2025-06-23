cmap = [
    0.000, 0.447, 0.741;  % blue
    0.929, 0.694, 0.125;  % yellow
    0.494, 0.184, 0.556;  % purple
    0.466, 0.674, 0.188;  % green
    0.301, 0.745, 0.933;  % cyan
    0.750, 0.550, 0.130;  % orange (less red than lines)
    0.200, 0.600, 0.500;  % teal
    0.635, 0.078, 0.941;  % violet/magenta (not red-toned)
];

set(groot, 'defaultAxesColorOrder', cmap)

clear cmap