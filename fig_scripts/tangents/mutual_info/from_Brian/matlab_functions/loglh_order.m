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
%
% modified 2/05 -- so that the input is a matrix of ISIs, where each row
% represents a sequence of n ISIs

% begin
clear
load fly199_isi

% create ISI distributions
for i = 1:12
    isi{i} = fly199_isi{i};
end
n_isi_dist = size(isi,2); % number of ISI vectors

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

% create prob density functions
for i = 1:n_isi_dist
    tmp = [];
    tmp = histc(log10(isi{i}),bins);
    tmp(find(tmp==0)) = 1; % ensure that each bin has at least one count;
    tmp2(i,:) = tmp./sum(tmp);
    clear tmp;
end

isi_dist{1} = tmp2(1:6,:);
isi_dist{2} = tmp2(7:12,:);
clear tmp2;

% maps elements of isi_mat to probability; the number of columns
% represents the number or repetitions.
for i = 1:size(isi_mat,1)
    for j = 1:size(isi_mat,2)
        % transform ISIs to bin number
        tmp = log10(isi_mat{i,j});
        tmp = floor((tmp-mn)./dx+2); % "+2" since "Inf's" are outside bins

        % bin numbers limited to the range of available bins
        % NOTE: this part of the code should never be needed...
        if any(find(tmp > (length(bins)-1))) || any(find(tmp <= 0))
            disp('Some data were outside of constructed bin vector');
        end
        tmp(find(tmp > (length(bins)-1))) = length(bins)-1;
        tmp(find(tmp <= 0)) = 1;

        % for prob(i,j) rows and columns correspond to indices for isi_mat
        prob{i,j}.px = reshape(isi_dist{j}(i,tmp),size(tmp,1), ...
                size(tmp,2));
        % use abs(j-3) to make index 1 when j = 2 and vice versa
        prob{i,j}.qx = reshape(isi_dist{abs(j-3)}(i,tmp),size(tmp,1), ...
                size(tmp,2));
        
        % row indicates Switching experiment number while column represents
        % whether data came from first distribution or second distribution
        D{i,j} = cumsum( log( (prob{i,j}.px) ./ ... 
                (prob{i,j}.qx) ));
        snr{i,j} = mean(D{i,j},2).^2 ./ ...
                var(D{i,j},0,2);
    end
end

%---------------------------------------------

for i = 1:6
    subplot(3,2,i)
    plot(snr{i,2},'r')
    hold on;
    plot(snr{i,1},'b')
    legend([num2str(var_key(i,1)) ' to ' num2str(var_key(i,2))], ...
        [num2str(var_key(i,2)) ' to ' num2str(var_key(i,1))])
end



