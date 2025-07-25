function [spike_times] = findspikes(V,thr,SampleRate)
% [spike_times] = findspikes(V,thr,SampleRate)
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
%
% NOTE: spikes are not counted if they occur within some time window such
% as 2 msec; this assumes that SampleRate is entered as Hz!
%
%-------------------------------------------------------------------------
% Brian Lundstrom, modified 0602, 0801

[M N] = size(V);
if N == 1
    V = V';
elseif N ~= 1 && M ~= 1
    error('Input V is not a 1 x N vector.')
end

%create row vector of spikes with 1's when V = threshold and 0's otherwise
%thr = 0; %set threshold for spike 'marker'
spikes = V >= thr;
spikes(1)=0; %set first number to zero in case spike is at beginning
spikes(length(spikes))=0; %set last number to zero in case spike is at end
spikes = diff(spikes);
%result is segments of 1 0 0 ... -1 where peak is somewhere in middle

%don't count any spikes within first 1 msec
spikes(1:round(SampleRate*0.001)) = 0;

start_peak = find(spikes == 1);
end_peak = find(spikes == -1); clear spikes;
if length(start_peak)>0
    if start_peak(1)>end_peak(1)
        end_peak = end_peak(2:end);
    end
end
if length(end_peak)==1 && length(start_peak)==0
    end_peak = [];
end
if length(start_peak) ~= length(end_peak)
    disp('Error in Spike Counting.')
end

spike_times = zeros(1,length(start_peak));

%for i=1:length(start_peak)
%    spike_times(i) = start_peak(i)/SampleRate;
%end

%don't count a spike if the mean of the previous msec > thr
%{
for i = 1:length(start_peak)
    if mean(V(start_peak(i)-round(SampleRate*0.001)+1:start_peak(i)))>thr-20
        start_peak(i) = -start_peak(i);
        end_peak(i) = -end_peak(i);
    end
end
start_peak = start_peak(start_peak>0);
end_peak = end_peak(end_peak>0);
%}

for i=1:length(start_peak)
    max_v = max(V(start_peak(i):end_peak(i)));
    tmp = find(V(start_peak(i):end_peak(i))==max_v);
    time_index = tmp(1); %take first value if >1 values exist, i.e. there are equal peak values
    %don't count a spike if the mean of the previous msec > thr
    %if mean(V(start_peak(i)-round(SampleRate*0.001)+1:start_peak(i)))<thr
        spike_times(i) = (start_peak(i)+time_index-1)/SampleRate;
    %end
end

if length(spike_times) >0
    %remove spike times that are within ABREF sec of the previous
    abref = 0.002;
    spike_times = spike_times([1 (find(diff(spike_times)>abref)+1)]);
    tmp = find(diff(spike_times)<abref);
    if isempty(tmp), tmp = 0;end
    disp([ num2str(tmp) ' spikes closer than ' num2str(1000*abref) ' msec found. Spikes closer than ' num2str(1000*abref) ' msec are not allowed.'])
end

return;