function [C_spike C_prior V D sproj sprior] = stc(Stim,Spikes,SampleRate,TimeSpan,Isolate)
%function [SpikeMatrix PriorMatrix V D sproj sprior] = stc(Stim,Spikes,SampleRate,TimeSpan,Isolate)
%
% Covariance analysis 
%
% INPUT
% Stim          =   stimulus (note: should not contain zero/initial point)
% Spikes        =   spike times (sec)
% SampleRate    =   sampling frequency of stimulus
% TimeSpan      =   range for stimulus window (sec), e.g. [-0.04 0.005]
% Isolate       =   0 = not isoloated; otherwise, use value in secs
%
% OUTPUT
% C_spike       =   covariance matrix of spike-conditioned stimulus
% C_prior       =   covariance matrix of all stimulus
% V             =   eigenvectors
% D             =   eigenvalues
% sproj         =   standardized projection onto stimuli (5 pos / 5 neg)
% sprior        =   standardized projection onto prior (5 pos / 5 neg)
%
% Covariance matrix:
% C = mean( spk_stim1*spk_stim2 ) - sta1*sta2 - mean( stim1*stim2 )
%
%----------------------------
% Brian Lundstrom, 2/05
% revised 4/05, 10/06

if nargin == 4; Isolate = 0; end;

%Use spikes that are preceeded by at least ISOLATE secs
tmp = [0 diff(Spikes)];
Spikes = Spikes(tmp > Isolate);

Stim = Stim./std(Stim);
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

C_spike = cov(SpikeMatrix);
C_prior = cov(PriorMatrix);
C_diff = C_spike - C_prior;
[V,D] = eig(C_diff);
[Y I] = sort(abs(diag(D)),'descend');

% hack to use 'positive' eigenvectors
for i  = 1:10
    if sum(V(:,I(i))) < 0
        V(:,I(i)) = -V(:,I(i));
    end
end

% scatter plot of stimulus projections along the first 5 eigenvectors
for i = [1:5 length(V)-4:length(V)]
    proj{i} = SpikeMatrix*V(:,i);
    prior{i} = PriorMatrix*V(:,i);
    std_prior(i) = std(prior{i});
    sproj{i} = proj{i}./std_prior(i);
    sprior{i} = prior{i}./std_prior(i);
end
