function plot_rate_isi(SpikeTimes, CycleTime, BinWidth)
% function plot_rate_isi(SpikeTimes, CycleTime, Title)
%
% INPUT
% SpikeTimes -- vector of spike times (sec)
% 
% BNL 2006-05
str = 'tmp';

% (1) Spike response for cycle
BinWidth = 0.1;
subplot(221)
[BinCenters,ResponsePerBin] = plot_spike_rate(SpikeTimes,CycleTime,... 
    BinWidth,0);
title(str,'FontSize',14)


% (2) Rate distributions
T = CycleTime;
Bins = min(ResponsePerBin):max(ResponsePerBin);
if length(Bins) <10, Bins = min(ResponsePerBin)-5:max(ResponsePerBin)+5; end;
P = CycleTime;
split = [0 P];
clear x y
for i = 1:length(split)-1
    a = split(i)/BinWidth + 1;
    z = split(i+1)/BinWidth;
    m = floor(T/BinWidth);
    [y(:,i) x(:,i)] = hist(ResponsePerBin(a:a+m-1),Bins);
    y(:,i) = y(:,i)/sum(y(:,i));
end

% Note:
% x = rate (spikes/sec)
% y = number of bins (from above) at the given rate

subplot(223)
set(gca,'FontSize',12)
plot(x,y,'.-')
xlabel('Rate (spikes/sec)')
ylabel('Rate density')
legend('\sigma_1','\sigma_2','\sigma_3','\sigma_2')
legend boxoff
title(['Data from first ' num2str(T) ' seconds of each epoch'])


% (4) ISI distributions
spikes = rem(SpikeTimes,CycleTime);
Bins = -3:0.5:1;
clear x2 y2
for i= 1:length(split)-1
    a = split(i);
    z = split(i+1);
    tmp = diff(spikes(find(spikes >= (a) & spikes <= (a+T))));
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi{i} = tmp(tmp2);
    [y2(:,i) x2(:,i)] = hist(log10(isi{i}),Bins);
    y2(:,i) = y2(:,i)/sum(y2(:,i));
end

subplot(224)
set(gca,'FontSize',12)
plot(x2,y2,'.-')
xlabel('Log_{10} \Delta (sec)')
ylabel('Interval density')
legend('\sigma_1','\sigma_2','\sigma_3','\sigma_2')
legend('Location','NorthWest')
legend boxoff
title(['Data from first ' num2str(T) ' seconds of each epoch'])

