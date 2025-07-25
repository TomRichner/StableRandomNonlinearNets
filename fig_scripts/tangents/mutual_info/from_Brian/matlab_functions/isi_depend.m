function mi = isi_depend(isi,dx,mn)

% Test independence of interpsike intervals (ISIs) by comparing the joint
% probability distribution of ISI_1 and ISI_n with the product of their
% separate probabilities. This comparison can be made using mutual
% information (see equ 4.11 in Dayan and Abott).
%
% In general, this code will create two matrices. The first matrix will
% consist of the joint probability density, where rows will designate the
% first ISI and columns will designate the second ISI. The second matrix
% will be similar expect that the elements will be the simple product of
% the individual probabilities.
%
%-------------------------------------------------------------------------
% Brian Lundstrom, 2/05

% load data
%load frb217_isi;
%isi = frb217_isi;

% begin
n_isi_dist = size(isi,2); % number of ISI vectors

% create bin vector
mx = 2;
%mn = -3;
%dx = 0.05;
edges = mn:dx:mx;
edges = [-inf edges inf];
tic;

% create prob density functions for single ISIs
for i = 1:1 %n_isi_dist
    tic;
    tmp = [];
    % BIN vector gives the bin number for each value of ISI{i}
    [tmp bin] = histc(log10(isi{i}),edges);
    % ensure that each "bin" of tmp has at least one count;
    tmp(find(tmp==0)) = 1; 
    isi_dist(i,:) = tmp./sum(tmp);
    clear tmp;

    % create joint prob density matrix, where row is first ISI and column is
    % nth ISI
    %nisi = 2; %for example, n=2 means look for two successive ISIs
    for nisi = 1:10
        
    % JOINT PROBABILITY = PDM1
    pdm1 = zeros(length(edges),length(edges));
    for m = 1:length(edges)
        for n = 1:length(edges)
            % Find a given bin value and then look N places after it in the
            % ISI vector
            tmp = find(bin(1:(length(bin)-(nisi-1))) == m) + (nisi-1);
            a = sum(bin(tmp) == n); %number of occurrences
            b = length(bin)-(nisi-1); %total possible
            pdm1(m,n) = a/b; %prob of the sequence
        end
    end

     
    % create prob density matrix of product of the two probabilities
    % NOTE: this step is equivalent to multiplying isi_dist'*isi_dist
    pdm2 = isi_dist(i,:)'*isi_dist(i,:);
%  
%    pdm2 = zeros(length(edges),length(edges));
%    for m = 1:length(edges)
%        for n = 1:length(edges)
%            a1 = sum(length(find(bin == m))); %number of occurrences
%            a2 = sum(length(find(bin == n))); %number of occurrences
%            b = length(bin);
%            pdm2(m,n) = a1*a2 / b^2; %product of probabilities
%        end
%    end

    %mutual information
    tmp = pdm1.*log2(pdm1./pdm2);
    tmp(isnan(tmp))=0; %replace NaN's with zeros
    mi(nisi,i) = sum(sum(tmp));
    end
    
    disp([num2str(i/n_isi_dist)*100 '% complete...'])
    toc;
end
toc;