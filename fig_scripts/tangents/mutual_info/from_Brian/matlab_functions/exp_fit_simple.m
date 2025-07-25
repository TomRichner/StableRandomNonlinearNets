function [A tau Fadapt sseBin Rf] = exp_fit(x,y,Period,Flag)
%function [A tau Fadapt sseBin ] = exp_fit(x,y,Period,Flag)
%
% Fits data to the equation: A * exp( t / tau),
%
% where Fadapt is:
%
% [initial firing rate - steady state FR] / inital FR
%
% Flag = 0 [no plot] or 1 [plot]
%-------------------------------------------------------------------------
% 
% BNL 2006-05

if nargin < 4
    Flag = 0;
end

if Flag ~=0
    figure
    %subplot(2,2,[1:2])
    plot(x,y), grid on, hold on
    xlabel('Time (sec)'); ylabel('Firing Rate (Hz)');
    hold on
end

% fit curves with exponential
xdata{1} = x(x<=Period/2);
xdata{2} = x(x>Period/2);
ydata{1} = y(x<=Period/2);
ydata{2} = y(x>Period/2);
y_offset(1) = mean(y(x>=Period/2-Period/11 & x<=Period/2)); %10 is last 4 bins of 15 per half period; 11-15 is last 3 bins
y_offset(2) = mean(y(x>=Period-Period/11 & x<=Period)); %7 is last 5 bins, 30 is 2 bins, 31 is 1
x_offset(1) = min(xdata{1});
x_offset(2) = min(xdata{2});

cp = { [1 0 0] [0 .5 0] };
for i = 1:2
    yy = ydata{i}-y_offset(i);
    xx = xdata{i}-x_offset(i);
    
    [estimates{i}, model{i}] = fitcurve_exp(xdata{i}-x_offset(i), ydata{i}-y_offset(i));
    %[estimates{i}, model{i}] = fitcurve_exp2(xdata{i}-x_offset(i), ydata{i});
    [sse(i) FittedCurve{i}] = model{i}(estimates{i});
    %subplot(2,2,i+2)
    if Flag ~=0
        %plot(xdata{i}, ydata{i}), hold on
        plot(xdata{i}, FittedCurve{i}+y_offset(i), 'Color', cp{i}), grid on
        %plot(xdata{i}, FittedCurve{i}, 'Color', cp{i}), grid on
    end
    Rf{i} = FittedCurve{i}+y_offset(i);
    A(i) = estimates{i}(1);
    tau(i) = (1/estimates{i}(2)); %time constant in sec
    Fadapt(i) = A(i)/y_offset(i)*100;
    sseBin(i) = sse(i)/length(xdata{i});
    tt{i} = char({['A = ' sprintf('%0.2g', A(i)) ' Hz (' sprintf('%0.3g',Fadapt(i)) '%),']; ...
        ['\tau = ' sprintf('%0.2g',tau(i)) ' sec, sse/bin = ' sprintf('%0.2g',sse(i)/length(xdata{i}))]});
end
if Flag ~=0
    h = legend(['Total bins = ' num2str(length(x))], tt{1},tt{2},'Location','SouthEast');
    legend boxoff
    set(h,'FontSize',12)
end