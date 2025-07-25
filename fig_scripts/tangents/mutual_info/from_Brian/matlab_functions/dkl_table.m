function dkl = dkl_table(px)
% Dkl = dkl_table(px);
%
% Calculation of Kullback-Leibler Divergence (Dkl).
%
% NOTE: Dkl is equal to calculation of D1 log-likelihood when n equals the
% population
%
% INPUT
% PX = probability distributions as column vectors
%
% OUTPUT
% DKL = Kullback-Leibler Divergence, where rows and columns represent the
% possible ISI distributions as:
%
% Dkl(r,c) = sum(px{r}.*log(px{r}./px{c}));
%
% Infinities and NaNs that result from the fraction log px(r)/px(c) will be
% ignored in the summation.
%
%-----------------------
% Brian Lundstrom, 8/04
% modified 12/04, 9/05

%load frb217_isi_2to10
%isi = frb217_isi;

s = size(px,2);

for r = 1:s
    for c = 1:s
        tmp = log2(px(:,r)./px(:,c));
        
        % ignore infinities and NaNs if present
        ind = isfinite(tmp);
        if length(ind)<length(tmp)
            disp(['Infinities and/or NaNs were ignored:' ...
            ' DKL results may be inaccurate.'])
        end
        
        dkl(r,c) = sum(px(ind,r).*tmp(ind));
    end
end
