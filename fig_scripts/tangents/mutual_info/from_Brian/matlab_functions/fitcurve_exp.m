function [estimates, model] = fitcurve_exp(xdata, ydata)
%function [estimates, model] = fitcurve_exp(xdata, ydata)
%
% code from MatLab help

% Call fminsearch with a random starting point.
start_point = [ydata(1)-ydata(end) 1/xdata(end)];
%start_point = [.1 1/xdata(end)];
model = @expfun;

estimates = fminsearch(model, start_point, ...
    optimset('MaxFunEvals',10000,'MaxIter',10000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve ErrorVector] = expfun(params)
        A = params(1);
        lambda = params(2);
        FittedCurve = A .* exp(-lambda * xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end