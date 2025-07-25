function iso_SpikeTimes = isolate_spikes(SpikeTimes, Tsilence)
%function iso_SpikeTimes = isolate_spikes(SpikeTimes, Tsilence)
%
% INPUTS
% SpikeTimes    =   vector of spiking times (sec)
% Tsilence      =   amount of time to separate spikes (sec)
%
% Throws out any spike times where the preceeding spike is less than
% T_silence before the current spike.
%
%------------------------------------------------
% BNL 05/05

loj_diff = logical([0 (diff(SpikeTimes) > Tsilence)]);
iso_SpikeTimes = SpikeTimes(loj_diff);