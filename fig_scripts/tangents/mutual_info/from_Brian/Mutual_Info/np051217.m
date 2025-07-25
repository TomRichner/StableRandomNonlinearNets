% variance estimation, etc
% run hh_var_est.m first

for ind = [1 14];

    std_true{ind} = E(round(est_t{ind}*SampleRate));
    stdu{ind} = unique(est_std{ind}); %unique values of estimated std

    for i = 1:length(stdu{ind})
        tmp = find(est_std{ind}==stdu{ind}(i));
        x{ind}(i) = mean(std_true{ind}(tmp));
        x_err{ind}(i) = std(std_true{ind}(tmp)); %/sqrt(length(std_true{ind}(tmp)));
    end
    
    errorbar(stdu{ind},x{ind},x_err{ind})
    hold on
end


% for rate
clear x x_err
for ind = [1 14];

    std_true{ind} = E(round(rt{ind}*SampleRate));
    stdu{ind} = unique(rate_std{ind});

    for i = 1:length(stdu{ind})
        tmp = find(rate_std{ind}==stdu{ind}(i));
        x{ind}(i) = mean(std_true{ind}(tmp));
        x_err{ind}(i) = std(std_true{ind}(tmp)); %/sqrt(length(std_true{ind}(tmp)));
    end
    
    errorbar(stdu{ind},x{ind},x_err{ind})
    hold on
end

% calculate correlations vs time window
for ind = 1:24
    std_true_int{ind} = E(round(est_t{ind}*SampleRate));
    std_true_rate{ind} = E(round(rt{ind}*SampleRate));
    int_corr(ind) = corr(std_true_int{ind}',est_std{ind}');
    rate_corr(ind) = corr(std_true_rate{ind}',rate_std{ind}');
end
plot(mean_win,int_corr,mean_win,rate_corr)

% calculate error bars for SNR calculations
sub = 11;
tic
for i = 1:30
    [snr D1 D2 pdm1 pdm2] = joint_ind_snr(sub);
    SNR{i} = snr;
    disp(i)
    toc
end

for i = 1:30
    [x y] = plot_snr_depend(SNR{i});
    joint_y(:,i) = y(:,1);
    ind_y(:,i) = y(:,2);
end
joint_err = std(joint_y,0,2)/sqrt(size(joint_y,2));
ind_err = std(ind_y,0,2)/sqrt(size(ind_y,2));


%-------------------------------------------------------------------------
% reconstruct D probability distributions for SNR plots
% data is from frb217_postn1000.mat
[prob D snr] = loglh_tmp(frb217_isi,data);

% adapted from np051021.m
h2l_ind = [D{1,3}(5,:) D{1,6}(5,:) D{1,9}(5,:) D{5,2}(5,:) D{5,7}(5,:) D{4,8}(5,:)];
l2h_ind = [D{3,1}(5,:) D{6,1}(5,:) D{9,1}(5,:) D{2,5}(5,:) D{7,5}(5,:) D{8,4}(5,:)];

%N = 50; % number of bins
x = -10:0.2:30;

figure
set(gca,'FontSize',12)
subplot(211)
tmp = hist(l2h_ind,x);
y_ind = tmp' / sum(tmp);
bar(x,y_ind)
title({'Frb217: considering 2 intervals, independent distributions'; ...
    'Low to High STD Change by 10x'})
subplot(212)
tmp = hist(h2l_ind,x);
y_ind(:,2) = tmp' / sum(tmp);
bar(x,y_ind(:,2))
title('High to Low STD Change by 10x')
