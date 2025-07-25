function SpikesParse = parse_spikes(Spikes, CycleTime, TotalTime)
%function SpikesParse = parse_spikes(Spikes, CycleTime, TotalTime)
%
% Split spike time vector into 2 parts according to a fixed cycle time such
% as during a 2-condition experiment
%
% INPUT:
% Spikes        =   spike time vector (sec)
% CycleTime     =   length of each condition (sec)
% TotalTime     =   time of trial, approx equal to max(Spikes)
%
% OUTPUT:
% SpikesParse{1}    =   spikes in first condition
% SpikesParse{2}    =   spikes in second condition
%
%------------------------------------------------
% BNL, 4/05

SpikesParse = cell(1,2);
delta = 2*CycleTime;
t = 0:delta:TotalTime;
for i = 1:length(t)
    tmp = Spikes( find( Spikes>t(i) & Spikes< (t(i)+CycleTime)) );
	tmp = tmp - t(i) + (i-1)*CycleTime;
    SpikesParse{1} = [SpikesParse{1} tmp];
end

t = CycleTime:delta:TotalTime;
for i = 1:length(t)
    tmp = Spikes( find( Spikes>t(i) & Spikes< (t(i)+CycleTime)) );
    tmp = tmp - t(i) + (i-1)*CycleTime;
    SpikesParse{2} = [SpikesParse{2} tmp];
end

