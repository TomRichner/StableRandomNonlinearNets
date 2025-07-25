close all
clear
clc


T = [0 100];
fs = 200;
dt = 1/fs;
t = T(1):dt:T(2);
nt = length(t);

[bH, aH] = butter(1,0.5/(fs/2),'high');

n_ch_in = 10;
rnd_signal = cumsum(randn(nt,n_ch_in));

u_ex = filtfilt(bH,aH,rnd_signal);

y = u_ex;

figure(1)
plot(u_ex)

%%

n_shift = 100;
mi_vec = zeros(n_shift,1);

chan_perm = randperm(n_ch_in)'
eq = chan_perm == (1:n_ch_in)'
% chan_perm = 1:n_ch_in;

y_flat = reshape(u_ex(:,chan_perm)',[],1);

u_ex_flat = reshape(u_ex',[],1);

for i_shift = 1:n_shift


    u_ex_flat_shifted = circshift(u_ex_flat,(i_shift-1)*n_ch_in,1);

    mi_vec(i_shift,1) = mi_fcn(u_ex_flat_shifted, y_flat);


end

figure(2)
plot(mi_vec)


% 
% mi_fcn(u_ex)
% 
% 
% n_shift = 100;
% mi_vec = zeros(n_shift,1);
% for i_shift = 1:n_shift
% 
%     B_shift=circshift(B,i_shift-1,1);
% 
%     mi_vec(i_shift,1) = mi_fcn(A,B_shift);
% 
% end
% 
% figure;
% plot(mi_vec)


