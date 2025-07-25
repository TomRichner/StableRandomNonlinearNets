
% multiple variances, zero mean
clear
vv = [400 ]; %900 1600 2500 3600];
for i = 1:length(vv);
        total_time = 10800;
        sample_rate = 20000; %hz
        stim_mean = 0;
        stim_var = vv(i);
        tic; 
        [SpikeTimes StimStats Stims v] = hhrk4_loop(total_time,stim_mean,stim_var,sample_rate);
        disp(['Repetition ' int2str(i)]); 
        fname = ['hh_' num2str(vv(i)) 'var'];
        save(fname,'SpikeTimes','StimStats','Stims')
        clear Stims StimStats SpikeTimes
        toc;
end
    
% multiple variances, multiple means
clear
vv = [0 100 400];
dc = [40 70];
for i = 1:length(vv);
    for j = 1:length(dc);
        total_time = 900;
        sample_rate = 20000; %hz
        stim_mean = dc(j);
        stim_var = vv(i);
        tic; 
        [SpikeTimes StimStats Stims Vtrace] = hhrk4_loop(total_time,stim_mean,stim_var,sample_rate);
        disp(['Repetition ' int2str(i)]); 
        fname = ['hh2_' num2str(dc(j)) 'dc_' num2str(vv(i)) 'var'];
        save(fname,'SpikeTimes','StimStats','Stims','Vtrace')
        clear Stims StimStats SpikeTimes Vtrace
        toc;
    end
end