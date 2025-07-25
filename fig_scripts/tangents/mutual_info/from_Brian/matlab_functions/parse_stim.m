function StimParse = parse_stim(Stim, TotalTime, ParseVector, SampleRate)
%function WavesParse = parse_stim(Stim, TotalTime, ParseVector, SampleRate)
%
% Split spike time vector into multiple parts according to a fixed cycle 
% time such as during a 2-condition experiment
%
% INPUT:
% Stim         =   stimulus time vector (sec)
% TotalTime    =   time of trial, approx equal to max(Spikes)
% SampleRate   =   sample rate (hz)
% ParseVector  =   signifies how to split stimulus vector,
%                   e.g. [0 8 12 20 24] for Period = 24 means to split the
%                   into 8, 4, 8, and 4 sec blocks                
%
% OUTPUT:
% StimParse{1}    =   stim in first condition
% StimParse{2}    =   stim in second condition, etc.
%
%------------------------------------------------
% BNL, 4/05, 0603

N = length(ParseVector)-1;
StimParse = cell(1,2);

Period = ParseVector(end);
row = Period*SampleRate;
col = TotalTime/Period;
tmp = reshape(Stim,row,col);

for i = 1:N
    StimParse{i} = reshape(tmp(ParseVector(i)*SampleRate+1: ...
        ParseVector(i+1)*SampleRate,:), ...
        diff(ParseVector(i:i+1))*SampleRate*col,1);
end
