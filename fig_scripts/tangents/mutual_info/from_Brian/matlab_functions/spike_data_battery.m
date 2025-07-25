function spike_data_battery(Stim, Spikes, Vtrace, SampleRate)
%function spike_data_battery(Stim, Spikes, Vtrace, SampleRate)
%
% This function yields plots for basic consistency checks of spike data
%
% INPUT
% Stim          =   stimulus wave form as cell array
% Spikes        =   spike times of response (sec) as cell array
% Vtrace        =   voltage trace of response as cell array
% SampleRate    =   sample rate (Hz) 
%
% OUTPUT PLOTS
% (1) time vs spikes
% (2) snippets of stim
% (3) snippets of vtrace
%
%------------------------------------------------
% BNL, 4/05

figure
color = 'brgcmykbrgcmyk';
for i = 1:length(Spikes)
    plot(Spikes{i},color(i))
    hold on;
    xlabel('nth spike')
    ylabel('Time')
    grid on;
end
legend('1','2','3','4')

for n = 1:length(Stim)
figure
index = length(Stim{n})/4;
for i = 1:4
    subplot(2,2,i)
    vec = floor(((i-1)*index+1)):(floor((i-1)*index+SampleRate*2));
    t = vec./SampleRate;
    plot(t,Stim{n}(vec))
    title(['Stim ' num2str(n)]);
    xlabel('seconds')
end
end

for n = 1:length(Vtrace)
    if isempty(Vtrace{n}); break; else
figure
index = length(Vtrace{n})/4;
for i = 1:4
    subplot(2,2,i)
    vec = floor(((i-1)*index+1)):(floor((i-1)*index+SampleRate*2));
    t = vec./SampleRate;
    plot(t,Vtrace{n}(vec))
    title(['Stim ' num2str(n)]);
    xlabel('seconds')
end
    end
end
