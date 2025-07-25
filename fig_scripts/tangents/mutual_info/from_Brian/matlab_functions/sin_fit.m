function [A phi sseBin ] = sin_fit(x,y,T,Flag,Tphi)
%function [A phi sseBin ] = sin_fit(x,y,T,Flag,Tphi)
%
% Fits data to the equation:
%       A*sin(2*pi/T*xdata+Tphi+phi)+mean(ydata)
%
% Flag = 0 if no plot is desired
%
%-------------------------------------------------------------------------
% 
% BNL 2007-06

if nargin < 4
    Flag = 0;
end

% fit curves with sine
% use last 

cp = { [1 0 0] [0 .5 0] };

[estimates, model] = fitcurve_sine(x, y, T, Tphi);
[sse FittedCurve] = model(estimates);

A = estimates(1); % Hz
phi = estimates(2); %phase (radians)
sseBin = sse/length(x);

if Flag ==1
    figure
    plot(x,y), grid on, hold on
    xlabel('Time (sec)'); ylabel('Firing Rate (Hz)');
    hold on

    plot(x, FittedCurve, 'Color', cp{2}), grid on

    tt = char({['A = ' sprintf('%0.2g', A) ' Hz (A/mu = ' sprintf('%0.2g', A/mean(y)) ') ' ];  ...
        ['\phi = ' sprintf('%0.2g',phi) ' Rad ']; ['sse/bin = ' sprintf('%0.2g',sse/length(x))]});
    h = legend(['Total bins = ' num2str(length(x))], tt,'Location','NorthEast');
    legend boxoff
    set(h,'FontSize',12)
end