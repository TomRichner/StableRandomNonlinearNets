function WavesParse = parse_waves(Waves, CycleTime, TotalTime, SampleRate)
%function WavesParse = parse_waves(Waves, CycleTime, TotalTime, SampleRate)
%
% Split spike time vector into 2 parts according to a fixed cycle time such
% as during a 2-condition experiment
%
% INPUT:
% Waves         =   stimulus time vector (sec)
% CycleTime     =   length of each condition (sec)
% TotalTime     =   time of trial, approx equal to max(Spikes)
% SampleRate    =   sample rate (hz)
%
% OUTPUT:
% SpikesParse{1}    =   spikes in first condition
% SpikesParse{2}    =   spikes in second condition
%
%------------------------------------------------
% BNL, 4/05

WavesParse = cell(1,2);
row = CycleTime*SampleRate;
col = TotalTime/CycleTime;
tmp = reshape(Waves,row,col);

WavesParse{1} = reshape(tmp(:,1:2:col),row*col/2,1);
WavesParse{2} = reshape(tmp(:,2:2:col),row*col/2,1);
