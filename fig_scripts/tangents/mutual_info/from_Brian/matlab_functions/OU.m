function Stim = OU(t, dt, mu, sigma, tau);
%Stim = OU(t, dt, mu, sigma, tau);
%
%Ornstein-Uhlenbeck for Giugliano paper
% t = total time
% dt = time step
% mu = mean of noise
% sigma = SD of noise
% tau = time constant of exponential filter

rand('state',sum(100*clock));
%randn('state',0);
Stim = zeros(1,ceil(t/dt));
R = randn(1,ceil(t/dt));

Stim(1) = mu;
for i = 2:ceil(t/dt)
    %These 2 forms are apparently equivalent -- the (1-exp(-x)) factor is
    %some sort of approx that gives the terms a slightly bigger value
    
    %Stim(i) = Stim(i-1) + (1 - exp(-dt/tau)) * (mu - Stim(i-1)) + sqrt( 1-exp(-2*dt/tau) ) * sigma * R(i-1);
    Stim(i) = Stim(i-1) +dt/tau * (mu - Stim(i-1)) +  sqrt(2*dt/tau) * sigma * R(i-1);
end