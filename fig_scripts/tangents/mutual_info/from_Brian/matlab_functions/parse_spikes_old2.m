function SpikesParse = parse_spikes(Spikes, Period, ParseVector)
%SpikesParse = parse_spikes(Spikes, Period, ParseVector)
%
% Split spike time vector into 2 parts according to a fixed cycle time such
% as during a 2-condition experiment
%
% INPUT:
% Spikes        =   spike time vector (sec)
% Period        =   length of both conditions (sec)
% ParseVector   =   signifies which part of each block to draw spikes from,
%                   e.g. [0 1; 5 6] means to include spikes from the first
%                   second and the last second of each condition
%
% OUTPUT:
% SpikesParse{:,1}    =   spikes in first condition
% SpikesParse{:,2}    =   spikes in second condition
%
%------------------------------------------------
% BNL, 0602

[M N] = size(ParseVector);
SpikesParse = cell(M,2);

for j = 1:2
    tmp = rem(Spikes,Period);
    offset = (j-1)*Period/2;
    for i = 1:M
        SpikesParse{i,j} = Spikes(find(tmp > ParseVector(i,1) + offset ...
            & tmp < ParseVector(i,2) + offset));
    end
end

% make a continuous spike time vector without jumps for the other condition

for j = 1:2
    for i = 1:M
        offset = (j-1)*Period/2;
        SpikesParse{i,j} = SpikesParse{i,j} - offset;
        SpikesParse{i,j} = floor(SpikesParse{i,j}/Period)*Period/2 + ...
            rem(SpikesParse{i,j},Period);
    end
end
