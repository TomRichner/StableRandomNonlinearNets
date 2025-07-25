
clear

total_time = 10;
spikes = [];
stats = [];
sample_rate = 20000; %hz
v = 2.3:0.1:4.1;
vv = [0 100 400 900];
dc = 0:5:500;

for i = 1:length(dc);
    for j = 1:4
        stim_mean = dc(i);
        tic; 
        stim_var = vv(i); %10^v(i);
        [spike_times stim_stat stim v] = hhrk4_loop(total_time,stim_mean,stim_var,sample_rate);
        Spikes{i} = spike_times;
        Stats{i} = stim_stat;
        Stims{i} = stim;
        Vtrace{i} = v;
        disp(['Repetition ' int2str(i)]); 
        toc;  
    end
end

save hh_data Spikes Stats
    