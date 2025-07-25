function [dkl] = isi_to_dkl(isi)
%function [dkl] = isi_to_dkl(isi)
%
% This function uses (my) standard values to create a table of the
% Kullback-Leibler divergence between distributions created from a cell
% array of ISI vectors.
%
% INPUT
% ISI       -   cell array of ISIs
%
% OUTPUT
% DKL       -   defined as Dkl(r,c) = sum(px{r}.*log(px{r}./px{c})); 
%
%-------------------------------------------------------------------------
% BNL, 05-9

mn = -4;
del = 0.05;
mx = 2;

[px nx bins] = isi_to_px(isi,mn,del,mx,'log10',1);
dkl = dkl_table(px);