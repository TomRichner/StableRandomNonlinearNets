% INT CODE MS
% 10/05
%
% based on np050920.m

% DKL comparison of distributions as a function of time from the variance
% switch. The steady-state comparison distribution is taked to be the last
% second before the switch for each variance condition.

% ISI distributions-------------------------------------------------------
% load data
str = 'hst23303'; %'rdb15901'
CycleTime = 24; %sec
BinWidth = 0.01; %sec
TimeStep = 1; %sec time step for x-axis of plot
xWidth = [1];
figure

for m = 1:length(xWidth)

SpikeTimes = read_bin(str);
P = CycleTime/6;
split = [0 2*P 3*P 5*P 6*P];
spikes = rem(SpikeTimes,CycleTime);

% make steady-state distributions for each of the 4 epochs
TP = 1; %seconds -- width of steady-state bin
for i = 1:length(split)-1
    z = split(i+1);
    tmp = diff(spikes(find(spikes > z-TP & spikes < z-TP+xWidth(m))));
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi_ss{i} = tmp(tmp2);
end

% make isi cell array for each epoch -- widening epoch from switch
%ct = 1;
%for i = 1:length(split)-1
%    a = split(i);
%    z = split(i+1);
%    steps = (z-a)/TimeStep;
%    for j = 1:steps
%        tmp = diff(spikes(find(spikes > a & spikes < a+j*TimeStep)));
%        % get rid of negative numbers at boundaries
%        tmp2 = logical((sign(tmp)+1) / 2);
%        isi{ct} = tmp(tmp2);
%        ct = ct+1;
%    end
%end

x = 0:TimeStep:CycleTime;
for i = 1:length(x)-1
    tmp = diff(spikes(find(spikes > x(i) & spikes < x(i)+xWidth(m)))); 
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi{i} = tmp(tmp2);
end

e{1} = [isi_ss{1}; isi'];
%tmp = [isi_ss{2}; isi_ss{4}];
%e{2} = [tmp; isi'];
e{2} = [isi_ss{2}; isi'];
e{3} = [isi_ss{3}; isi'];
e{4} = [isi_ss{4}; isi'];

for i = 1:length(e)
    px = isi_to_px(e{i},-3,0.05,0,'log10',1);
    tmp = dkl_table(px);
    dkl{m}(:,i) = tmp(2:end,1); %1 since isi_ss are in first cell of e{i}
end

subplot(2,2,m+2)
t = 0:TimeStep:CycleTime-TimeStep;
set(gca,'FontSize',12)
plot(t,dkl{m},'.-')
legend('SS \sigma_1','SS \sigma_2','SS \sigma_3','SS \sigma_2')
title({['ISI DKL for increasing bins after switch: ' str] ; ... 
    [' with xWidth = ' num2str(xWidth(m)) ' sec']})
xlabel('Seconds')
grid on
legend boxoff
end


% Rate distributions------------------------------------------------------
% load data
str = 'hst23303'; % 'rdb15901';
CycleTime = 24; %sec
BinWidth = 0.01; %sec
SpikeTimes = read_bin(str);
TimeStep = .3; %sec time step for x-axis of plot
xWidth = [.3];

for m = 1:length(xWidth)

P = CycleTime/6;
split = [0 2*P 3*P 5*P 6*P];
spikes = rem(SpikeTimes,CycleTime);

% bin spiketime data in order to find spike rates
EndTime = floor(max(SpikeTimes)/CycleTime)*CycleTime;
SpikeTimes = SpikeTimes(1:max(find(SpikeTimes<EndTime)));
[fr BinCenters] = hist(spikes, CycleTime/BinWidth);
ResponsePerBin = fr.*1/BinWidth.*CycleTime./EndTime;
Bins = floor(min(ResponsePerBin)):ceil(max(ResponsePerBin));

% Steady-state is last 2 seconds of epoch
TP = 1; %seconds -- width of steady-state bin
for i = 1:length(split)-1
    z = split(i+1)/BinWidth;
    tmp = hist(ResponsePerBin(z-TP/BinWidth+1: ... 
        z-TP/BinWidth+xWidth(m)/BinWidth),Bins);
    % ensure that each bin has at least one count;
    tmp(find(tmp==0)) = 1;
    rate_ss{m}(:,i) = tmp/sum(tmp);
end

x = 0:TimeStep:CycleTime;
for i = 1:length(x)-1
    tmp = hist(ResponsePerBin(round(x(i)/BinWidth)+1: ... 
        round(x(i)/BinWidth+xWidth(m)/BinWidth)),Bins);
    % ensure that each bin has at least one count;
    tmp(find(tmp==0)) = 1;
    rate{m}(:,i) = tmp/sum(tmp);
end

for i = 1:size(rate_ss{m},2)
    for j = 1:size(rate{m},2)
        dkl_rate{m}(j,i) = sum(rate{m}(:,j).*log(rate{m}(:,j)./rate_ss{m}(:,i)));
    end
end

t = 0:TimeStep:CycleTime-TimeStep;
set(gca,'FontSize',12)
subplot(2,2,m)
plot(t,dkl_rate{m},'.-')
legend('SS \sigma_1','SS \sigma_2','SS \sigma_3','SS \sigma_2')
title({['Rate DKL for increasing bins after switch: ' str] ; ... 
    [' with xWidth = ' num2str(xWidth(m)) ' sec']})
xlabel('Seconds')
grid on
legend boxoff

end

