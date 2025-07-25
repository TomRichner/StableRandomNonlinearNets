% Generate HH ISI distributions

clear

TotalTime = 1800;
LoopTime = 180;
StimMean = 0;
SampleRate = 20000;
StimVar = [ 10993 17428 ]; % [276 438 694 1100 1744 2761 4375 6938 10993 17428];
%StimVar = [276 357 438 566 694 897 1100 1422 1744 2253 2761 3568 4375 5656 6938 8965 10993 14210 17428 20000];


for i = 1:length(StimVar)
    [SpikeTimes StimStats Stims Vtrace] = hhrk4_loop(TotalTime,LoopTime,StimMean,StimVar(i),SampleRate);
    rate(i) = length(SpikeTimes)/TotalTime;
    stim_var(i) = mean(StimStats(:,2));
    hh_isi{i} = diff(SpikeTimes);
    i/length(StimVar)
end

