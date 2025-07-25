function [mu eb r n] = circmean(phi, alpha)

% mu = circ_mean(phases)
%   Computes the mean direction for circular data.
%
%   Input:
%     phi	sample of angles in radians in columns
%     alpha for confidence interval, e.g. 0.05--> 95% CI
%
%   Output:
%     mu	mean direction
%     eb    error bar widtch, upper/lower confidence limit
%     r     resultant vector lenght for each mu
%     n     number of samples
%
% References:
%   Biostatistical Analysis, J. H. Zar
%
% See also (parts of code adated from):
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
%
% B Lundstrom, Winter 2009

tmp = size(phi);
if tmp(1) == 1
    phi = phi';
end

% find sum of cos and sin of angles
rsum = nansum(exp(1i*phi));
n = sum(~isnan(phi));

% find mean and r
mu = angle(rsum);
r = abs(rsum)./n;

% find confidence limits, Zar 26.7
R = n.*r;
c2 = chi2inv((1-alpha),1);

% check for resultant vector length and select appropriate formula
for i = 1:length(r)
    if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
      tmp = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
    elseif r(i) >= .9
      tmp = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
    else 
      tmp = NaN;
      %warning('Resultant vector does not allow to specify confidence limits on mean. \nResults may be wrong or inaccurate.');
    end
    % apply final transform
    eb(i) = acos(tmp/R(i));
end
