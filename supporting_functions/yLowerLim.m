function out = yLowerLim(newLower)
% yLowerLim  Get or set the lower y-axis limit.
%   curr = yLowerLim()         returns the current lower y-limit.
%   yLowerLim(newLower)        sets the lower limit to newLower,
%                              keeping the upper limit unchanged.
%
% Example:
%   plot(randn(100,1));
%   yLowerLim(-2)              % force the bottom of the y-axis to −2

    yl = ylim;
    if nargin == 0
        out = yl(1);
    else
        ylim([newLower, yl(2)]);
    end
end
