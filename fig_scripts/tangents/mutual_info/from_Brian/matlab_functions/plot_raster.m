function plot_raster(SpikeTimes, CycleTime)
%function plot_raster(SpikeTimes, CycleTime)
%
% BNL, 05/09 

xdata = rem(SpikeTimes,CycleTime);
ydata = ceil(SpikeTimes./CycleTime);
scatter(xdata,ydata,'.');