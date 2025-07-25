function [STAy STAx NLy NLx Prior Posterior] = sta_threshold(Stim,Spikes,SampleRate,TimeSpan,P)
% [STAy STAx NLy NLx] = sta_threshold(Stim,Spikes,SampleRate,TimeSpan)
% Find Spike-triggered average and the nonlinear threshold function
%
% INPUT
% Stim          =   stimulus (note: should not contain zero/initial point)
% Spikes        =   spike times (sec)
% SampleRate    =   sampling frequency of stimulus
% TimeSpan      =   range for stimulus window (sec), e.g. [-0.04 0.005]
% P             =   0 for no plot, 1 for plot
%
%-------------------------------------------------------------------------
% Brian Lundstrom, 0602

Stim = Stim/std(Stim);
Stim = Stim - mean(Stim);

win = abs(TimeSpan) .* SampleRate; %indicates length of window prior to spike time
SpikeIndex = round(Spikes .* SampleRate); % map spike times to stimulus index
tmp = SpikeIndex > win(1) & SpikeIndex < (length(Stim)-win(2)+1); %don't use partial windows
SpikeIndex = SpikeIndex(find(tmp)); %spikes = indices of spike times >= win
PriorIndex = win(1)+1 + floor( rand(1,length(SpikeIndex)).*(length(Stim)-sum(win)-1) );
SpikeMatrix = zeros(length(SpikeIndex),sum(win)+1);
PriorMatrix = zeros(length(PriorIndex),sum(win)+1);

% create matrix of spike-triggered stimuli
for i = 1:length(SpikeIndex)
    SpikeMatrix(i,:) = Stim((SpikeIndex(i)-win(1)):(SpikeIndex(i)+win(2)))';
    PriorMatrix(i,:) = Stim((PriorIndex(i)-win(1)):(PriorIndex(i)+win(2)))';
end

STAy = mean(SpikeMatrix);

% threshold function P(spike|x1) = [ P(spike) P(x1|spike) ] / P(x1)
%   where P(spike) ~ rate, P(x1|spike) is dist. of conditional stimuli 
%   projected on f1/f2 and P(x1) is prior

prior_proj = (PriorMatrix*STAy')/std(PriorMatrix*STAy');
%post_proj = (SpikeMatrix*STAy')/std(SpikeMatrix*STAy');
post_proj = (SpikeMatrix*STAy')/std(PriorMatrix*STAy');
rg = max(abs(post_proj));
nbins = round(length(Spikes)/50);
%del = 2*rg/(nbins-1);
bins = linspace(-rg,rg,nbins);

% add ones to all bins -- equivalent to assuming Bayesian uniform prior
Prior = (hist(prior_proj,bins)+1)./(length(Spikes)+nbins);
Posterior = (hist(post_proj,bins)+1)./(length(Spikes)+nbins);
Likelihood = Posterior./Prior/sum(Posterior./Prior);

NLy = cumsum(Likelihood);
NLx = bins;
STAx = TimeSpan(1):1/SampleRate:TimeSpan(2);

if P==1
    figure
    subplot(2,1,1)
    plot(STAx,STAy)
    xlabel('Time (sec)')
    ylabel('Stimulus (SD)')
    title('Spike-triggered stimulus')
    subplot(2,1,2)
    plot(bins,NLy,bins,Prior,bins,Posterior)
    legend('Cum Prob(Spike|Stim)','Prob(Stim)','Prob(Stim|Spike)')
    legend('Location','NorthWest')
    legend boxoff
    title('Cumulative Prob(Spike|Stimulus)')
    xlabel('Stimulus (SD)')
    ylabel('Probability of firing')
end



