function [SpikeTimes StimStats Stims Vtrace] = hhrk4_loop(TotalTime,LoopTime,StimMean,StimVar,SampleRate)
%function [SpikeTimes StimStats Stims Vtrace] = hhrk4_loop(TotalTime,LoopTime,StimMean,StimVar,SampleRate)
%
% Hodgkin and Huxley model of neuronal membrane
%
% INPUT 
% TotalTime     =   total time of simulation in seconds
% LoopTime      =   time per loop -- loops are needed to save memory (600s
%                   max)
% StimMean      =   desired mean of stimulus
% StimVar       =   desired variance of stimulus
% SampleRate    =   sampling rate (Hz)
%
% OUTPUT
% SpikeTimes    =   1 x n vector of spike times
% StimStats     =   actual mean and variance of stimulus
% 
%-------------------------------------------------------------------------
% Brian Lundstrom, Autumn 2004
% revised 5/05


% Initial values:
    Vi = 0;
    ni = 0.3177;
    mi = 0.0529;
    hi = 0.5961;
    yi = [Vi ni mi hi];
    dt = 1000/SampleRate; %msec
    SpikeTimes = [];
    Vtrace = [];
    stimbit1 = [];
    v_tot = [];
    Stims = [];
    down_sample = SampleRate/4000; %downsample to specificed Hz

    %LoopTime = 600; %seconds of data per loop
    tspan = [0 LoopTime*1000];
    loopn = ceil(TotalTime/LoopTime);
    disp([num2str(loopn*LoopTime/60) ' minutes of data will be created.'])

% loop
for i = 1:loopn
    
% Stimulus    
    stim_length = (LoopTime+1/SampleRate)*2;
    [stim stimbit2] = make_stim(stim_length,StimMean,StimVar,SampleRate,stimbit1);
    stimbit1 = stimbit2;
    StimStats(i,:) = [mean(stim) var(stim)];
    tmp = stim(1:(2*down_sample):length(stim)); %downsample by factor
    Stims = [Stims tmp(2:length(tmp))]; %remove initial point
    clear tmp;

% Solve HH equations via RK4 solver    
    [tt yy] = hhrk4(stim, tspan, dt, yi);
    if max(yy(1,:))>200 || min(yy(1,:))<-200 || sum(isnan(Vtrace))>0
        disp('************************Voltage potential is out of bounds.')
    end

% Find time for spikes
    st = findspikes(yy(1,:),0,20000) + (i-1)*LoopTime;
    SpikeTimes = [SpikeTimes st];
    
    tmp = yy(1,1:down_sample:length(yy(1,:)));
    Vtrace = [Vtrace tmp(2:length(tmp))];
    yi = yy(:,length(yy))';
	clear tt yy stim tmp;
    
end