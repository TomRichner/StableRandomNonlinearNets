function padYLim(x)
% padYLim  Extend current y-axis limits by ±x of the current range.
%   padYLim(x) adds fractional padding x (e.g. 0.1 for 10%) to both ends.
%
% Example:
%   plot(randn(100,1));
%   padYLim(0.2);   % adds 20% padding to top and bottom

    ylim(ylim + [-x x]*diff(ylim));
end
