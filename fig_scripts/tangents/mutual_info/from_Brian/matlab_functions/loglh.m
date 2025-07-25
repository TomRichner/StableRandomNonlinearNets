function [prob D snr Dm] = loglh(isi)
% [prob D snr] = loglh(isi)
%
% Calculates the log likelihood ratio using 2 vectors of ISIs.
%
% Log likelihood ratio
%   = sum( log( P(isi_1/distrA)/P(isi_1/distrB) + log( P(isi_2/distrA ...
%
% EXPLANATION:
% A series of ISIs from the 1st data set is compared to the constructed
% probability density function from each of the data sets to determine the
% respective probability of that given ISI occurring.
%
% NOTE:
% Resulting values are dependent on two practical factors: 
%   (1) number of repetitions or times that the ISI data is sampled, i.e. 
%   REP_MAX or the number of D's calculated for a given n
%   (2) the size of the histogram bin size, since the routine ensures that
%   all bins have at least one ISI (otherwise log 0 = inf)
%
% plotting -- see NP040812.M
%
% INPUT
% ISI = an array of ISI vectors
%
% OUTPUT
% REPS = structure array containing fields:
%   DATA = ISI indices
%   PROB = probabilites corresponding to the ISI indices
%   D = log-likelihood discriminations
%   SNR = calculated as mean(D)^2 / var(D)
%
% For the arrays PROB, D, and SNR, the row indicates where the data came
% from while the column value indicates which distribution is being used
% for the comparison.
%
%-----------------------
% Brian Lundstrom, 8/04 v2
% modified 12/04
% modified 09/05

% begin
n_isi_dist = size(isi,2); % number of ISI vectors
n = 100; % max number of ISIs considered
reps = 1000

% create probability mass functions
mx = 2;
mn = -3;
dx = 0.05;
[isi_dist nx bins] = isi_to_px(isi,mn,dx,mx,'log10',11);

% create data vectors and map indices of data vector to bin number and then
% to probability; each column is a sequence of n numbers, where the first
% number is a randomly picked 'seed' ISI. Therefore, the number of columns
% represents the number or repetitions.
for i = 1:n_isi_dist
    
    % construct table of ISI indices for a given distribution/ISI vector
    tmp = floor((length(isi{i})-n-1)*rand(1,reps))+1;
    data{i} = ones(n,1)*tmp + [0:(n-1)]'*ones(1,reps); 
    
	% transform ISIs to bin number
    tmp = log10(isi{i}(data{i}));
    tmp = floor((tmp-mn)./dx+2); % "+2" since "Inf's" are outside bins
    
    % bin numbers limited to the range of available bins
    % NOTE: this part of the code should never be needed...
    if any(find(tmp > (length(bins)-1))) || any(find(tmp <= 0))
        disp('Some data were outside of constructed bin vector');
    end
    tmp(find(tmp > (length(bins)-1))) = length(bins)-1;
    tmp(find(tmp <= 0)) = 1;
    
    % for prob{i,j} row indicates where data came from, column indicates 
    % comparison distribution, i.e. P( Del | Dist) == P ( i | j)
    %
    % in a given prob{i,j} array, rows are the interval number and columns
    % are the sample number
    for j = 1:n_isi_dist
        prob{i,j} = reshape(isi_dist(tmp,j),n, reps);
    end
    
    % row indicates P(x) and column is Q(x)
    for j = 1:n_isi_dist
        D{i,j} = cumsum( log( prob{i,i} ./ prob{i,j} ),1);
        snr{i,j} = mean(D{i,j},2).^2 ./ var(D{i,j},0,2);
    end
end

% the expectation of D1 is the Kullback-Leibler divergence
for i = 1:size(snr,1)
    for j = 1:size(snr,2)
        for n = 1:size(D{i,j},1)
        Dm{n}(i,j) = mean(D{i,j}(n,:));
        end
        d_var(i,j) = var(D{i,j}(1,:));
        snr1(i,j) = snr{i,j}(5);
    end
end
Dm
snr1
d_var
%

%---------------------------------------------

% plotting -- see NP040812.M


