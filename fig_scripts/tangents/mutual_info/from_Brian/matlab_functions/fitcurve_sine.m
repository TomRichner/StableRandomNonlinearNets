function [estimates, model] = fitcurve_sine(xdata, ydata,T,Tphi)
%function [estimates, model] = fitcurve_sine(xdata, ydata)
%
% T = period of sine wave
%
% adapted from code from MatLab help (fitcurve_exp)

% Call fminsearch with a random starting point.
start_point = [20 -.90];
model = @sinfun;

estimates = fminsearch(model, start_point, ...
    optimset('MaxFunEvals',10000,'MaxIter',10000));
% sinfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = sinfun(params)
        A = params(1);
        phi = params(2);
        FittedCurve = A*sin(2*pi/T*xdata+phi+Tphi)+mean(ydata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end