function [estimates, model] = fitcurve_exp2(xdata, ydata)
%function [estimates, model] = fitcurve_exp(xdata, ydata)
%
% code from MatLab help

% Call fminsearch with a random starting point.
guess(1) = ydata(1)-ydata(end);
guess(2) = 1/(xdata(1)-xdata(end));
guess(3) = mean(ydata);

start_point = guess;
model = @expfun;

estimates = fminsearch(model, start_point, ...
    optimset('MaxFunEvals',10000,'MaxIter',10000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun(params)
        A = params(1);
        lambda = params(2);
        B = params(3);
        FittedCurve = A .* exp(-lambda * xdata) + B;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end