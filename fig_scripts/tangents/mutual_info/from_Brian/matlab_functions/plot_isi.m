function plot_isi(isi)
% function plot_isi(isi)
%
% Plots an array of ISIs as a histogram, such that the sum of all (x,y)
% pairs equals 1. Bins span the range log10[-4 2] seconds.
%
% BNL

if ~iscell(isi)
    isi{1} = isi;
end

Bins = -4:0.02:2;
for i = 1:length(isi)
    [y(:,i)] = hist(log10(isi{i}),Bins);
    y(:,i) = y(:,i)/sum(y(:,i));
end

plot(Bins,y)
xlabel('Log_{10} \Delta (sec)')
ylabel('Interval density')