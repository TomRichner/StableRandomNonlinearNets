function [estimates, model] = fitcurve_hill(xdata, ydata, base)
%function [estimates, model] = fitcurve_sine(xdata, ydata)

% Call fminsearch with a random starting point.
Starting = [max(ydata)-min(ydata) .5*(max(ydata)-min(ydata))+min(ydata) 1.1 min(ydata)];
model = @hillfun;

estimates = fminsearch(model, Starting,...
    optimset('MaxFunEvals',10000,'MaxIter',10000));
% hillfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = hillfun(params)
       %Hill equation: base + Rmax x^h / (R50^h + x^h)
        Rmax=params(1);
        D50=params(2);
        h=params(3);
        base = params(4);
        FittedCurve = base + Rmax*xdata.^h./ (D50^h + xdata.^h);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end