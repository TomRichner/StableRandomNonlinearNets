%-------------------------------------------------------------------------
% HHRK4 -- Plot ISI distributions from ISI vector array
%
clear
load hhrk4_isi;
isi = hhrk4_isi;
n_isi_dist = size(isi,2); % number of ISI vectors

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

% create prob density functions
for i = 1:n_isi_dist
    tmp(i,:) = histc(log10(isi{i}),bins);
    isi_dist(i,:) = tmp(i,:)/sum(tmp(i,:));
    clear tmp;
end

% create ISI distributions for hhrk4 data
color = 'bgrcmykbgrcmyk';
mat = 1:n_isi_dist;
for i = 1:length(mat)
    plot(bins(10:120),isi_dist(mat(i),10:120),color(i))
    hold on;
end
xlabel('Log10 of ISI Interval')
ylabel('Probability Density Normalized Over Bin Total, Bin Width = 0.05s')
title('HHRK4 ISI Distributions')
legend('1','2','3','4','5','6','7','8','9','10')

figure;
for i = 1:length(mat)
    plot(10.^(bins(10:65)),log10(isi_dist(mat(i),10:65)),color(i))
    hold on;
end
xlabel('ISI Interval')
ylabel('Log10 of Probability Density Normalized Over Bin Total, Bin Width = 0.05s')
title('HHRK4 ISI Distributions')
legend('1','2','3','4','5','6','7','8','9','10')

%-------------------------------------------------------------------------
% FRB217 -- Plot ISI distributions from ISI vector array
%
clear
load frb217_isi_2to10;
isi = frb217_isi;
n_isi_dist = size(isi,2); % number of ISI vectors

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

% create prob density functions
for i = 1:n_isi_dist
    tmp(i,:) = histc(log10(isi{i}),bins);
    isi_dist(i,:) = tmp(i,:)/sum(tmp(i,:));
    clear tmp;
end

% create ISI distributions for hhrk4 data
color = 'bgrcmykbgrcmyk';
mat = 1:n_isi_dist;
for i = 1:length(mat)
    plot(bins(10:120),isi_dist(mat(i),10:120),color(i))
    hold on;
end
xlabel('Log10 of ISI Interval')
ylabel('Probability Density Normalized Over Bin Total, Bin Width = 0.05s')
title('FRB217 ISI Distributions')
legend('Experiment 2','3','4','5','6','7','8','9','10')
