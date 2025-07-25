function Spikes_bin = isi_to_spikes(isi,del)
%function Spikes_bin = isi_to_spikes(isi,del)
%
% Create a 0's and 1's vector from ISI vector for a given delta T.
%
% INPUT:
% ISI - interspike-interval matrices
% DEL - width of each bin
%
%-------------------------------------------------------------------------
% BNL 0604

if ~iscell(isi)
    tmp = isi;
    clear isi;
    isi{1} = tmp;
end

for i = 1:length(isi)
    tmp = cumsum(isi{i});
    mx = max(tmp);
    tmp = ceil(tmp/del);
    Spikes_bin{i} = zeros(ceil(mx/del),1);
    Spikes_bin{i}(tmp) = 1;
end
