function [snr D1 D2 pdm1 pdm2] = joint_ind_snr(sub)

% create SNR based on joint probability

% begin
%clear
load frb217_isi_2to10
isi = frb217_isi;

% create ISI distributions
n_isi_dist = size(isi,2); % number of ISI vectors
rep = 1000;

% create isi probability distributions
mx = 0;
mn = -3;
dx = 0.05;
%sub = 0; % what to use to ensure non-zero bins
edges = mn:dx:mx;
edges = [-inf edges inf];
[px nx BinEdges BinLabels] = isi_to_px(isi,mn,dx,mx,'log10',sub);
tic

% create joint and independent probability distributions
for i = 1:n_isi_dist

    % create joint prob density matrix, where row is first ISI and column is
    % nth ISI
    %nisi = 2; %for example, n=2 means look for two successive ISIs
    nisi = 2;
    
    % BIN vector gives the bin number for each value of ISI{i}
    bin = BinLabels{i};
    pdm1{i} = zeros(length(BinEdges),length(BinEdges));
    for m = 1:length(BinEdges)
        for n = 1:length(BinEdges)
            tmp = find(bin(1:(length(bin)-(nisi-1))) == m) + (nisi-1);
            a = sum(bin(tmp) == n); %number of occurrences
            b = length(bin)-(nisi-1); %total possible
            pdm1{i}(m,n) = (a+1)/b; %prob of the sequence
        end
    end
    pdm1{i} = pdm1{i}./sum(sum(pdm1{i}));
    % create prob density matrix of product of the two probabilities
    % NOTE: this step is equivalent to multiplying isi_dist'*isi_dist
    pdm2{i} = px(:,i)*px(:,i)';
end

% create data vectors and map indices of data vector to bin number and then
% to probability; each column is a sequence of n numbers, where the first
% number is a randomly picked 'seed' ISI. Therefore, the number of columns
% represents the number or repetitions.
for i = 1:n_isi_dist
    
    % construct table of ISI indices for a given distribution/ISI vector
    data{i} = floor((length(isi{i})-nisi-1)*rand(1,rep))+1;
    data{i} = ones(nisi,1)*data{i} + ...
        [0:(nisi-1)]'*ones(1,rep); 
    
	% transform ISIs to bin number
    tmp = log10(isi{i}(data{i}));
    tmp = floor((tmp-mn)./dx+2); % "+2" since "Inf's" are outside bins
    
    % bin numbers limited to the range of available bins
    % NOTE: this part of the code should never be needed...
    if any(find(tmp > (length(edges)-1))) || any(find(tmp <= 0))
        disp('Some data were outside of constructed bin vector');
    end
    tmp(find(tmp > (length(edges)-1))) = length(edges)-1;
    tmp(find(tmp <= 0)) = 1;
    
    % for prob(i,j) row indicates where data came from, column indicates 
    % comparison distribution
    % prob1 = joint distribution
    % prob2 = independent distributions
    for j = 1:n_isi_dist
        prob1{i,j} = diag(pdm1{j}(tmp(1,:),tmp(2,:)))';
        prob2{i,j} = diag(pdm2{j}(tmp(1,:),tmp(2,:)))';
    end
    
    % row indicates P(x) and column is Q(x)
    for j = 1:n_isi_dist
        D1{i,j} = log( (prob1{i,i}) ./ ... 
            (prob1{i,j} ) );
        
        %D1{i,j}(isinf(D1{i,j}))=0;
        %tmp = D1{i,j}(find(D1{i,j}));
        
        snr_joint(i,j) = mean(D1{i,j},2).^2 ./ ...
            var(D1{i,j},0,2);
        
        D2{i,j} = log( (prob2{i,i}) ./ ... 
            (prob2{i,j} ) );
        
        %D2{i,j}(isinf(D2{i,j}))=0;
        %tmp = D2{i,j}(find(D2{i,j}));
          
        snr_ind(i,j) = mean(D2{i,j},2).^2 ./ ...
            var(D2{i,j},0,2);
    end
end
toc;
snr{1} = snr_joint;
snr{2} = snr_ind;
