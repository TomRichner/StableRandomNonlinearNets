function [SpikeMatrix PriorMatrix C_diff V D] = stca(Stim,Spikes,SampleRate,TimeSpan)
%function [SpikeMatrix PriorMatrix C_diff V D] = stca(Stim,Spikes,SampleRate,TimeSpan)
%
% Covariance analysis 
%
% INPUT
% Stim          =   stimulus (note: should not contain zero/initial point)
% Spikes        =   spike times (sec)
% SampleRate    =   sampling frequency of stimulus
% TimeSpan      =   range for stimulus window (sec), e.g. [-0.04 0.005]
%
% OUTPUT
% C_spike       =   covariance matrix of spike-conditioned stimulus
% C_prior       =   covariance matrix of all stimulus
% C_diff        =   C_spike - C_diff
%
% Covariance matrix:
% C = mean( spk_stim1*spk_stim2 ) - sta1*sta2 - mean( stim1*stim2 )
%
% NOTE: Stim and Spikes vectors should be one continuous condition since
% the prior distribution is taken from the Stim vector.
%
%----------------------------
% Brian Lundstrom, 2/05
% revised 4/05

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

C_spike = cov(SpikeMatrix);
C_prior = cov(PriorMatrix);
C_diff = C_spike - C_prior;
[V,D] = eig(C_diff);
[Y I] = sort(abs(diag(D)),'descend');

% hack to use 'positive' eigenvectors
for i  = 1:3
    if sum(V(:,I(i))) < 0
        V(:,I(i)) = -V(:,I(i));
    end
end


% plot covariance matrices
figure
t = TimeSpan(1):1/SampleRate:TimeSpan(2);
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
plot(t,V(:,I(1)),t,V(:,I(2)),t,V(:,I(3)))
title('Eigenvectors (sorted by magnitude)')
legend('1','2','3','Location','Northwest')
xlabel('Time (sec)')


% scatter plot of stimulus projections along the first 2 eigenvectors
for i = 1:3
    proj{i} = SpikeMatrix*V(:,I(i));
    prior{i} = PriorMatrix*V(:,I(i));
    std_prior(i) = std(prior{i});
    sproj{i} = proj{i}./std_prior(i);
    sprior{i} = prior{i}./std_prior(i);
end

%{
figure
set(gca,'FontSize',14)
scatter3(sprior{1},sprior{2},sprior{3},'.');
hold on;
scatter3(sproj{1},sproj{2},sproj{3},'r.');
xlabel('Projection on Eigenvector 1')
ylabel('Projection on Eigenvector 2')
zlabel('Projection on Eigenvector 3')
%}

figure
for i = 1:2
subplot(3,1,i)
scatter(sprior{1},sprior{i+1},'.')
hold on;
scatter(sproj{1},sproj{i+1},'r.')
xlabel('EV 1 Proj (SD Prior)')
ylabel(['EV ' num2str(i+1) ' Proj (SD Prior)'])
end

%{
subplot(3,1,3)
scatter(sprior{2},sprior{3},'.')
hold on;
scatter(sproj{2},sproj{3},'r.')
xlabel('EV 2 Proj (SD Prior)')
ylabel('EV 3 Proj (SD Prior)')
%}

%
% decision function P(spike|x1) = [ P(spike) P(x1|spike) ] / P(x1)
%   where P(spike) ~ rate, P(x1|spike) is dist. of conditional stimuli 
%   projected on f1/f2 and P(x1) is prior

bins = -4:.1:4;
rate = length(Spikes)/Spikes(length(Spikes)); % Hz
for i = 1:2
    pt{i} = hist(sproj{i},bins)./length(Spikes);
    p{i} = hist(sprior{i},bins)./length(Spikes);
    df{i} = pt{i}./p{i};
    tmp = isnan(df{i});
    df{i}(tmp) = 0;
end

subplot(3,1,3)
plot(bins,df{1},bins,df{2})
title('Prob(Spike|Stimulus) for 0.1 bin width')
legend('Projection onto eigenvector 1','Projection onto eigenvector 2')
xlabel('Stimulus')
ylabel('Probability Distribution of firing')
%}
