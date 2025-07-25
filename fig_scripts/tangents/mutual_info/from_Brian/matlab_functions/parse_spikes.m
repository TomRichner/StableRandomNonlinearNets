function SpikesParse = parse_spikes(Spikes, Period, ParseVector)
%SpikesParse = parse_spikes(Spikes, Period, ParseVector)
%
% Split spike time vector into multiple parts according to a fixed cycle 
% time such as during a 2-condition experiment
%
% INPUT:
% Spikes        =   spike time vector (sec)
% Period        =   length of n non-repeated conditions (sec)
% ParseVector   =   signifies which part of each block to draw spikes from,
%                   e.g. [0 1; 5 6; 6 7; 11 12] means to include spikes 
%                   from the first second and the last second of each 
%                   condition for T = 12 sec with 2 conditions
%
% OUTPUT:
% SpikesParse{:,1}    =   spikes in first condition
% SpikesParse{:,2}    =   spikes in second condition, etc.
%
%------------------------------------------------
% BNL, 0602, revised 0603

[M N] = size(ParseVector);
SpikesParse = cell(1,M);

tmp = rem(Spikes,Period);
for i = 1:M
    SpikesParse{i} = Spikes((tmp > ParseVector(i,1)...
        & tmp < ParseVector(i,2)));
end

% make a continuous spike time vector without jumps for other conditions
for i = 1:M
    SpikesParse{i} = floor(SpikesParse{i}/Period)*diff(ParseVector(i,:)) ...
        + rem(SpikesParse{i},Period) - ParseVector(i,1);
end
