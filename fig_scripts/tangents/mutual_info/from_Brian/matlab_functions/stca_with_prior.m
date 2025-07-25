function stca(Stim,Spikes,SampleRate,Range,StimPrior)
%
% Covariance analysis 
%
% Stim          =   stimulus
% Spikes        =   spike times (sec)
% SampleRate    =   sampling frequency of stimulus
% Range         =   range for stimulus window (sec)
% StimPrior     =   stimulus for prior -- used when spikes are parsed
%
% Covariance matrix:
% C = mean( spk_stim1*spk_stim2 ) - sta1*sta2 - mean( stim1*stim2 )
%
%----------------------------
% Brian Lundstrom, 2/05
% revised 4/05

%clear
%Spikes = st;
%Stim = s./std(s);
%SampleRate = 20000; % Hz sample rate
%Range = [-0.030 0.005]; % msec range for stimulus window

win = abs(Range) .* SampleRate; %indicates length of window prior to spike time

SpikeIndex = round(Spikes .* SampleRate); % map spike times to stimulus index
tmp = SpikeIndex > win(1) & SpikeIndex < (length(Stim)-win(2)); %don't use partial windows
SpikeIndex = SpikeIndex(find(tmp)); %spikes = indices of spike times >= win
PriorIndex = win(1)+1 + floor(rand(1,length(SpikeIndex)) .* (length(Stim)-sum(win)-1));
SpikeMatrix = zeros(length(SpikeIndex),sum(win));
PriorMatrix = zeros(length(PriorIndex),sum(win));

% create matrix of spike-triggered stimuli
for i = 1:length(SpikeIndex)
    SpikeMatrix(i,:) = Stim((SpikeIndex(i)-win(1)+1):(SpikeIndex(i)+win(2)))';
    PriorMatrix(i,:) = StimPrior((PriorIndex(i)-win(1)+1):(PriorIndex(i)+win(2)))';
end

C_spike = cov(SpikeMatrix);
C_prior = cov(PriorMatrix);
C_diff = C_spike - C_prior;
[V,D] = eig(C_diff);

% plot covariance matrices
figure
t = (Range(1)+1/SampleRate):1/SampleRate:Range(2);
subplot(3,2,1)
pcolor(t,t,C_spike); shading('flat');
title('Spike-triggered Covariance Matrix')
xlabel('Time (sec)')
ylabel('Time (sec)')
subplot(3,2,2)
pcolor(t,t,C_prior); shading('flat');
title('Prior Covariance Matrix')
subplot(3,2,3)
pcolor(t,t,C_diff); shading('flat');
title('Covariance Difference Matrix')
subplot(3,2,4)
plot(t,mean(SpikeMatrix))
title('Spike-triggered Average')
ylabel('Stimulus (Standard Deviations)')
subplot(3,2,5)
plot(diag(D),'.');
title('Eigenvalues of Covariance Difference Matrix')
subplot(3,2,6)
plot(t,V(:,1),t,V(:,2))
title('First 2 Eigenvectors')
legend('Eigenvector 1',' Eigenvector 2','Location','Northwest')
xlabel('Time (msec)')


% scatter plot of stimulus projections along the first 2 eigenvectors
figure
proj1 = SpikeMatrix*V(:,1);
proj2 = SpikeMatrix*V(:,2);
prior1 = PriorMatrix*V(:,1);
prior2 = PriorMatrix*V(:,2);
subplot(3,1,1)
scatter(proj1,proj2,'.')
title('2-D Distribution of Conditional Stimuli')
subplot(3,1,2)
scatter(prior1,prior2,'.')
title('2-D Distribution of Prior Stimuli')
xlabel('Projection onto Eigenvector 1')
ylabel('Projection onto Eigenvector 2')


% decision function P(spike|x1) = [ P(spike) P(x1|spike) ] / P(x1)
%   where P(spike) ~ rate, P(x1|spike) is dist. of conditional stimuli 
%   projected on f1/f2 and P(x1) is prior

%rate = 5/500; % prob of spike for each 2 msec bin
nbin = 30; % approximate number of bins
tmp = [min(proj1) min(proj2) max(proj1) max(proj2)];
bins = floor(min(tmp)):(max(tmp)-min(tmp))/nbin:floor(max(tmp));
t1 = hist(proj1,bins);
t1 = t1/sum(t1);
p1 = hist(prior1,bins);
p1 = p1/sum(p1);
t2 = hist(proj2,bins);
t2 = t2/sum(t2);
p2 = hist(prior2,bins);
p2 = p2/sum(p2);
df1 =  t1 ./ p1;
df2 =  t2 ./ p2;

set(gca,'FontSize',12)
subplot(3,1,3)
plot(bins,df1,bins,df2)
title('Decision functions for low variance condition')
legend('Projection onto eigenvector 1','Projection onto eigenvector 2')
xlabel('Stimulus')
ylabel('Probability of firing')


