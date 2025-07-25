T = [0 100];
fs = 1000;
dt = 1/fs;
t = T(1):dt:T(2);
nt = length(t);

[bH, aH] = butter(3,1/(fs/2),'high');

rnd_signal = cumsum(randn(nt,1));

u_ex = filtfilt(bH,aH,rnd_signal);

figure(1)
plot(u_ex)

n_shift = 100;
mi_vec = zeros(n_shift,1);
for i_shift = 1:n_shift

    u_ex_shifted = circshift(u_ex,i_shift-1,1);

    mi_vec(i_shift,1) = mi_fcn(u_ex_shifted, u_ex);

end

figure(10)
plot(mi_vec)

