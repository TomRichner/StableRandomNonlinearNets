function [MI] = mutual_info_SISO(data_in, data_out, n_bins_in, n_bins_out)
% computes mutual information. Deals with multiple inputs and multiple outputs by pooling their histograms
   % data_in and data_out shoudl be channels x time

    [n_ch_in, n_t_in]   = size(data_in);
    [n_ch_out, n_t_out] = size(data_out);
    
    assert(n_t_in == n_t_out, 'data_in and data_out must have same number of time samples');
    assert(n_ch_in == 1, 'data_in must have a single channel for SISO');
    assert(n_ch_out == 1, 'data_out must have a single channel for SISO');

    min_in = min(data_in,[],'all');
    max_in = max(data_in,[],'all');

    min_out = min(data_out,[],'all');
    max_out = max(data_out,[],'all');

    edges_in = linspace(min_in, max_in, n_bins_in + 1);
    edges_out = linspace(min_out, max_out, n_bins_out + 1);

    P_in = histcounts(data_in, edges_in);
    P_in = P_in ./ sum(P_in, 'all');

    P_out = histcounts(data_out, edges_out);
    P_out = P_out ./ sum(P_out, 'all');

    P_joint = histcounts2(data_in, data_out, edges_in, edges_out);
    P_joint = P_joint ./ sum(P_joint, 'all');

    MI = sum(P_joint .* log2(P_joint ./ (P_in' * P_out)), 'all', 'omitnan');

end