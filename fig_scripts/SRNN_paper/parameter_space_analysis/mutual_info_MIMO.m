function [MI] = mutual_info_MIMO(data_in, data_out, n_bins_in, n_bins_out, excluded_output_chs)
% computes mutual information. Deals with multiple inputs and multiple outputs by pooling their histograms
   % data_in and data_out shoudl be channels x time

    if not(isempty(excluded_output_chs))
        data_out(excluded_output_chs,:) = []; % exlude specific output channels
    end

    [n_ch_in, n_t] = size(data_in);
    [n_ch_out, n_t_out] = size(data_out);

    if n_t < n_ch_in
        error('n_t < n_ch_in, data_in should be channels by time')
    end
    if n_t < n_ch_out
        error('n_t < n_ch_out, data_out should be channels by time')
    end

    if n_t ~= n_t_out
        error('data in and data out time dimension does not match')
    end

    min_in = min(data_in,[],'all');
    max_in = max(data_in,[],'all');

    min_out = min(data_out,[],'all');
    max_out = max(data_out,[],'all');

    edges_in = linspace(min_in, max_in, n_bins_in+1);
    edges_out = linspace(min_out, max_out, n_bins_out+1);

    P_in_by_channel = zeros(n_ch_in, n_bins_in);
    for i_ch_in = 1:n_ch_in
        P_in_by_channel(i_ch_in,:) = histcounts(data_in(i_ch_in,:), edges_in);
    end
    P_in = sum(P_in_by_channel, 1) ./ sum(P_in_by_channel, 'all'); % sum across channels in and normalize to 1

    P_out_by_channel = zeros(n_ch_out, n_bins_out);
    for i_ch_out = 1:n_ch_out
        P_out_by_channel(i_ch_out,:) = histcounts(data_out(i_ch_out,:),edges_out);
    end
    P_out = sum(P_out_by_channel,1)./sum(P_out_by_channel,'all'); % sum across channels out and normalize to 1

    P_joint_sum = zeros(n_bins_in, n_bins_out);
    for i_ch_in = 1:n_ch_in
        for i_ch_out = 1:n_ch_out
            P_joint_sum = P_joint_sum + histcounts2(data_in(i_ch_in,:), data_out(i_ch_out,:), edges_in, edges_out);
        end
    end
    P_joint = P_joint_sum ./ sum(P_joint_sum, 'all');

    MI = sum(P_joint .* log2(P_joint ./ (P_in' * P_out)), 'all', 'omitnan');

end