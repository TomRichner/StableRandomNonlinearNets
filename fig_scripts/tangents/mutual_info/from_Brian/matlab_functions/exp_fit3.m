function [A tau Fadapt] = exp_fit3(x,y,Period)
%function [A tau Fadapt sseBin ] = exp_fit(x,y,Period)
%
% Fits data to the equation: A * exp( t / tau),
%
% where Fadapt is:
%
% [initial firing rate - steady state FR] / inital FR
%
%-------------------------------------------------------------------------
% 
% BNL 2006-05

figure
%subplot(2,2,[1:2])
plot(x,y), grid on, hold on
xlabel('Time (sec)'); ylabel('Firing Rate (Hz)');
hold on

% fit curves with exponential
% use last 

xdata{1} = x(x<Period/2);
xdata{2} = x(x>Period/2);
ydata{1} = y(x<Period/2);
ydata{2} = y(x>Period/2);
y_offset(1) = mean(y(x>Period/2-Period/4 & x<Period/2));
y_offset(2) = mean(y(x>Period-Period/4 & x<Period));
x_offset(1) = min(xdata{1});
x_offset(2) = min(xdata{2});

cp = { [1 0 0] [0 .5 0] };
for i = 1:2
    %[estimates{i}, model{i}] = fitcurve_exp(xdata{i}-x_offset(i), ydata{i}-y_offset(i));
    [estimates{i}, model{i}] = fitcurve_exp3(xdata{i}-x_offset(i), ydata{i});
    [FittedCurve{i}] = model{i}(estimates{i},xdata{i}-x_offset(i));
    %subplot(2,2,i+2)
    %plot(xdata{i}, ydata{i}), hold on
    %plot(xdata{i}, FittedCurve{i}+y_offset(i), 'Color', cp{i}), grid on
    plot(xdata{i}, FittedCurve{i}, 'Color', cp{i}), grid on
    A(i) = estimates{i}(1);
    tau(i) = (1/estimates{i}(2)); %time constant in sec
    Fadapt(i) = A(i)/y_offset(i)*100;
    
    tt{i} = char({['A = ' sprintf('%0.2g', A(i)) ' Hz (' sprintf('%0.3g',Fadapt(i)) '%),']; ...
        ['\tau = ' sprintf('%0.2g',tau(i)) ]});
end

h = legend(['Total bins = ' num2str(length(x))], tt{1},tt{2},'Location','SouthEast');
legend boxoff
set(h,'FontSize',12)