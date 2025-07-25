function Dkl_isi = dkl_table(isi)
% Dkl = dkl_table(isi);
%
% Calculation of Kullback-Leibler Divergence (Dkl). Bins are created with a
% log10 range with 20 bins per power of 10, e.g. from 1 microsecond to 1000
% seconds.
%
% NOTE: changes in the bin definitions result in some change to Dkl values.
% NOTE: Dkl is equal to calculation of D1 log-likelihood when n equals the
% population
%
% INPUT
% ISI = cell array of individual ISI vectors
%
% OUTPUT
% DKL = Kuellback-Leibler Divergence, where rows and columns represent the
% possible ISI distributions as:
%
% Dkl_isi(r,c) = sum(px{r}.*log(px{r}./px{c}));
%
%-----------------------
% Brian Lundstrom, 8/04
% modified 12/04, 9/05

%load frb217_isi_2to10
%isi = frb217_isi;
s = length(isi);

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

for i = 1:s
    px{i} = histc(log10(isi{i}),bins);
    px{i}(find(px{i}==0)) = 1; % ensure that each bin has at least one count;
    px{i} = px{i} / sum(px{i}); % normalization
end

% NOTE: Dkl is computed here using log not log2
for r = 1:s
    for c = 1:s
        Dkl_isi(r,c) = sum(px{r}.*log(px{r}./px{c}));
    end
end
