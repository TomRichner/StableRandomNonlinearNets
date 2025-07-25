function [spike_times] = findspikes(V,thr, SampleRate)
% [spike_times] = findspikes(V,thr, SampleRate)
%
% INPUTS
% V = 1xn membrane potential vector
% thr = threshold for spike counting
% SampleRate (Hz)
%
% OUTPUT
% SPIKE_TIMES = times in seconds corresponding to membrane potential maxima
%
% This routine finds areas of membrane voltage above a given threshold THR
% and then returns the times of the maxima of these areas

% use 2 msec window to search for spike max
win = round(SampleRate*0.002);
spike_times = [];

i = 2;
while i < length(V)-win
    if V(i) > thr
        tmp = max(V(i:i+win));
        tmp = find(V(i:i+win) == tmp);
        tmp = i+tmp(1)-1;
        spike_times = [spike_times tmp/SampleRate];
        i = i+win;
    end
    i = i+1;
end

return;