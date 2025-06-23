function out = yUpperLim(newUpper)
% yUpperLim  Get or set the upper y-axis limit.
%   curr = yUpperLim()         returns the current upper y-limit.
%   yUpperLim(newUpper)        sets the upper limit to newUpper,
%                              keeping the lower limit unchanged.
%
% Example:
%   plot(randn(100,1));
%   yUpperLim(2)               % force the top of the y-axis to +2

    yl = ylim;
    if nargin == 0
        out = yl(2);
    else
        ylim([yl(1), newUpper]);
    end
end
