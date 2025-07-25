function [estimates, model] = fitcurve_exp3(xdata, ydata)
%function [estimates, model] = fitcurve_exp(xdata, ydata)
%
% code from MatLab help

% Call fminsearch with a random starting point.
guess(1) = ydata(1)-ydata(end);
guess(2) = 100; %1/(xdata(1)-xdata(end));
guess(3) = mean(ydata);

start_point = guess;
model = @expfun;
xdata = xdata';
ydata = ydata';

[estimates,r,J] = nlinfit(xdata,ydata,model,start_point);
 
    function [FittedCurve sse] = expfun(params,xdata)
        A = params(1);
        lambda = params(2);
        B = params(3);
        FittedCurve = A .* exp(-lambda * xdata) + B;
    end
end