function yLimPadMin(x, maxYLowerLim, minYUpperLim)
% padYLimBounded  Extend current y-axis limits by ±x fraction of the current range,
%                 with optional bounds on the lower and upper limits.
%
%   padYLimBounded(x)  
%       pads both ends of the y-axis by x*range.
%
%   padYLimBounded(x, maxYLowerLim)  
%       first enforces yl(1) ≤ maxYLowerLim, then pads by x*newRange.
%
%   padYLimBounded(x, maxYLowerLim, minYUpperLim)  
%       enforces both yl(1) ≤ maxYLowerLim and yl(2) ≥ minYUpperLim, 
%       then pads by x*newRange.
%
% Inputs:
%   x             – fractional padding (e.g. 0.1 for ±10%)
%   maxYLowerLim  – (optional) the maximum allowable lower limit
%   minYUpperLim  – (optional) the minimum allowable upper limit
%
% Example:
%   plot(randn(100,1));
%   padYLimBounded(0.2, -1, 2);
%   % ensures lower ≤ –1, upper ≥ 2, then adds 20% padding.

    % get current limits
    yl = ylim;

    % enforce lower bound if provided
    if nargin >= 2 && ~isempty(maxYLowerLim)
        yl(1) = min(yl(1), maxYLowerLim);
    end

    % enforce upper bound if provided
    if nargin >= 3 && ~isempty(minYUpperLim)
        yl(2) = max(yl(2), minYUpperLim);
    end

    % compute new range and pad
    r         = diff(yl);
    ylPadded  = yl + [-x, x] * r;

    % apply
    ylim(ylPadded);
end
