% variance estimation, etc
% run hh_var_est.m first


est_t{n},est_std{n},rt{n},rate_std{n}

%isi -- find bias, variance of estimated std
ind = [1 14];
for i = 1:length(ind);
    std_true{i} = E(round(est_t{ind(i)}*SampleRate));
    isi_means(1,i) = mean(std_true{i});
    isi_means(2,i) = mean(est_std{i});
end
x = [2 15];

% for rate
for i = 1:length(ind);
    std_true{i} = E(round(rt{ind(i)}*SampleRate));
    rt_means(1,i) = mean(std_true{i});
    rt_means(2,i) = mean(rate_std{i});
end

%------------------------------------------------------------
% variance estimation, etc
% run hh_var_est.m first

ind = [1:24];
for i = 1:length(ind);
    std_true{ind(i)} = E(round(est_t{ind(i)}*SampleRate));
    mn_std = min(std_true{ind(i)});
    mx_std = max(std_true{ind(i)});
    bins = mn_std: (mx_std-mn_std)/9: mx_std;
    for m = 1:length(std_true{ind(i)})
        [y yi] = min(abs(bins-std_true{ind(i)}(m)));
        std_true{ind(i)}(m) = bins(yi);
    end
    
    for j = 1:length(bins)
        tmp = find(std_true{ind(i)}==bins(j));
        isi_m(j,i) = mean(est_std{ind(i)}(tmp));
        isi_s(j,i) = std(est_std{ind(i)}(tmp)); %/sqrt(length(est_std{ind(i)}(tmp)));
    end  
    vvar(i,1) = mean((est_std{ind(i)}-std_true{ind(i)}).^2);
end
x = repmat(bins',1,length(ind));
errorbar(x,isi_m,isi_s);




ind = [1:24];
for i = 1:length(ind);
    std_true{ind(i)} = E(round(rt{ind(i)}*SampleRate));
    mn_std = min(std_true{ind(i)});
    mx_std = max(std_true{ind(i)});
    bins = mn_std: (mx_std-mn_std)/9: mx_std;
    for m = 1:length(std_true{ind(i)})
        [y yi] = min(abs(bins-std_true{ind(i)}(m)));
        std_true{ind(i)}(m) = bins(yi);
    end
    
    for j = 1:length(bins)
        tmp = find(std_true{ind(i)}==bins(j));
        rate_m(j,i) = mean(rate_std{ind(i)}(tmp));
        rate_s(j,i) = std(rate_std{ind(i)}(tmp)); %/sqrt(length(rate_std{ind(i)}(tmp)));
    end     
    vvar(i,2) = mean((rate_std{ind(i)}-std_true{ind(i)}).^2);
end
x = repmat(bins',1,length(ind));
errorbar(x,rate_m,rate_s);

isi_bias = sum(abs(isi_m-x));
rate_bias = sum(abs(rate_m-x));

