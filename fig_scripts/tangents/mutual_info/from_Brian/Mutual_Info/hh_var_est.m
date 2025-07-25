% Hodgkin and Huxley model of neuronal membrane
%
% Brian Lundstrom, Nov 2005

clear

%__________________________________________________________
%
% create probability distributions from previous steady state ISI data
load hh_isi
stim_var = [276 438 694 1100 1744 2761 4375 6938 10993 17428];
%load hh_isi20;
%stim_var = [276 357 438 566 694 897 1100 1422 1744 2253 2761 3568 4375 ... 
%    5656 6938 8965 10993 14210 17428 20000];
isi = hh_isi;
mn = -4;
del = 0.05;
mx = 2;
dist = 'log10';
sub = 1;
[px nx bins] = isi_to_px(isi,mn,del,mx,dist,sub);

%__________________________________________________________
%
% make stimulus

% parameters
tspan = [0 60e3]; % msec
dt = 0.05; % msec
steps = dt/2; % use dt/2 as step size since this is required by RK4 solver
time = tspan(1):steps:tspan(2);
SampleRate = 40000;

% make stimulus; STIMBIT2 is used when concatenating stimulus vectors
% NOTE: stimulus length needs to be sampled at 2x of expected output since
% RK4 solvers sample at the half time step

StimLength = tspan(2)/1000;
StimMean = 0;
StimVar = 1;
SampleRate = 40000;
Stim = make_stim(StimLength,StimMean,StimVar,SampleRate);

%{
rand('state',sum(100*clock));
(length(Stim)-1)/100;
ups = floor(4*SampleRate);
e = ceil(StimLength./4);
tmp = randn(1,e+1);
%tmp = [1.9164 -0.1053 -0.2167 0.2979 -1.4508 1.4492 0.5734 -1.6486 ... 
    %-0.3790 0.7585 -0.8167 0.1587 0.2015 0.3809];
E = resample(tmp,ups,1);
E = E/std(E);
%E = E*sqrt(max(stim_var))/max(abs(E)); %scale E by max stim std
E = 25*E+65;
%}

load slow_stim051123
NewStim = E(1:length(Stim)).*Stim;
ts = (1:length(NewStim))/SampleRate; %time in sec


%W = randn(1, stim_length);
%tau = 6; %sec
%ts = tau*1000;
%E = zeros(1,stim_length);

% dS/dt = -1/tau S(t) + W(t), where W(t) is white noise
% S = conv(exp(-s/t), W)
%for i = 2:stim_length
%    E(i) = E(i-1) + dt/2 * ( -1/ts * E(i-1)) + W(i-1)*sqrt(dt/ts)*4;
%end
%
%stim = 3*E.*W*20;
%t = dt/2:dt/2:length(stim)*dt/2;

% measure STD of stimulus every est_del
est_del = SampleRate/2;
ind = 1:est_del:length(NewStim)-est_del;
for i = 1:length(ind)
    ss(i) = std(NewStim(ind(i):ind(i)+est_del));
    tt(i) = (ind(i)+est_del) / SampleRate; %time in sec 
end

% use Hodgkin-Huxley Runge-Kutta 4th order solver with appropriate initial
% values
[t y] = hhrk4(NewStim,tspan,dt,[0 0.3177 0.0529 0.5961]);
[spike_times] = findspikes(y(1,:),65,20000);
isi = [ 0 diff(spike_times)];

clear rate_std
isi_win = [1:1:24]; %number of isis-1 to consider in variance estimation
for n = 1:length(isi_win)
    ind = isi_win(n)+2:10:length(isi);
    for i = 1:length(ind)
        aa = ind(i)-isi_win(n);
        zz = ind(i);
        clear est
        for m = 1:size(px,2)
            %map ISIs to a bin number for the probability distributions
            bin_num = floor((log10(isi(aa:zz))-mn)/del+2);
            est(m) = sum(log(px(bin_num,m)));
        end
        [tmp yi] = max(est);
        est_std{n}(i) = sqrt(stim_var(yi));
        est_t{n}(i) = spike_times(zz);
        est_rate{n}(i) = (isi_win(n)+2)/(spike_times(zz)-spike_times(aa));
        est_std_time{n}(i) = (spike_times(zz)-spike_times(aa));
    end

    %calculate rate with fixed window
    mean_win(n) = mean(est_std_time{n});
    rt{n} = mean_win(n):mean_win(n):tspan(2)/1000; %defines right edge of bin
    tmp = hist(spike_times,[0 rt{n}]) / mean_win(n);
    tmp = tmp(1:end-1);
    for i = 1:length(tmp)
        [y yi] = min(abs(rate-tmp(i)));
        rate_std{n}(i) = sqrt(stim_var(yi));
    end
end

%{
for n = 1:length(isi_win)
    rt{1} = rt{1}(1:2:end);
    rate_std{1} = rate_std{1}(1:2:end);

    subplot(length(isi_win),1,n)
    plot(tt,ss,est_t{n},est_std{n},rt{n},rate_std{n},'.')
    legend('Stim Std','ISI Est','Rate Est')
    title({[num2str(isi_win(n)+1) ' ISIs Considered']; ... 
        ['Mean time for estimation = ' num2str(mean_win(n)) ' sec ']})
end
%}

% see also NP051217.m