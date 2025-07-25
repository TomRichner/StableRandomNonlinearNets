% INT CODE MS FIGURE
% based on NP050907.m
%
% Oct 2005, BNL


% examining hst233, 050908
str = 'hst23303';
CycleTime = 24; %sec
BinWidth = 0.01; %sec
SpikeTimes = read_bin(str);

% (2) Spike response for cycle
figure(2)
subplot(2,2,1)
[BinCenters,ResponsePerBin] = plot_spike_rate(SpikeTimes,CycleTime,... 
    BinWidth);
title(str,'FontSize',14)


% (3) Rate distributions -- take last 1 second
T = 1;
Bins = min(ResponsePerBin):max(ResponsePerBin);
P = CycleTime/6;
split = [0 2*P 3*P 5*P 6*P];
clear x y
for i = 1:length(split)-1
    a = split(i)/BinWidth + 1;
    z = split(i+1)/BinWidth;
    [y(:,i) x(:,i)] = hist(ResponsePerBin(z-T/BinWidth+1:z),Bins);
    y(:,i) = y(:,i)/sum(y(:,i));
end

% Note:
% x = rate (spikes/sec)
% y = number of bins (from above) at the given rate

subplot(2,2,3)
set(gca,'FontSize',12)
plot(x,y,'.-')
xlabel('Rate (spikes/sec)')
ylabel('Rate density')
legend('\sigma_1','\sigma_2','\sigma_3','\sigma_2')
legend boxoff


% (4) ISI distributions -- take last 1 second
spikes = rem(SpikeTimes,CycleTime);
Bins = -3:0.05:0;
clear x2 y2 isi
for i= 1:length(split)-1
    a = split(i);
    z = split(i+1);
    tmp = diff(spikes(find(spikes > (z-T) & spikes < (z))));
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi{i} = tmp(tmp2);
    [y2(:,i) x2(:,i)] = hist(log10(isi{i}),Bins);
    y2(:,i) = y2(:,i)/sum(y2(:,i));
end

subplot(2,2,4)
set(gca,'FontSize',12)
plot(x2,y2,'.-')
xlabel('Log_{10} \Delta (sec)')
ylabel('Interval density')
legend('\sigma_1','\sigma_2','\sigma_3','\sigma_2')
legend boxoff

