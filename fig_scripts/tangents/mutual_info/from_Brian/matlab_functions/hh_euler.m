function [t y ie] = hh_euler(Stim, SampleRate)

% Hodgkin-Huxley model neuron
%
% INPUT:
%  Stim          vector of stimulus input; length of output will be determined 
%                by length of stimulus divided by sampling frequency; note 
%                that stimulus is automatically filtered such that noise
%                standard deviation is maintained (i.e. filter is
%                normalized by its norm)
%
%  SampleRate    in kHz; determines time step for Euler Method; this should be at
%                least 20 kHz; 50 kHz may be preferable
%
% OUTPUT:
%   T            1 x n vector for time points in msec
%   Y            n x 4 vector; column 1 is voltage; columns 2-4 are the n, m,
%                and h gating variables
%   ie           external current injected into neuron
%
% EXAMPLE:
%   SampleRate = 50;
%   Stim = [zeros(1,25e3) 3*ones(1,50e3) 4.6*ones(1,50e3) 70*(rand(1,50e3)-0.5)];
%   [t y ie] = hh_euler(Stim, SampleRate);
%
% Equations:
% (1)   im = Gl(V - El) + Gk n^4(V - Ek) + GNa m^3 h(V - ENa)
% (2)   cm * (dV/dt) = -im + Ie/A
% (3)   dn/dt = alphan*(1-n) - betan*(n)
% (4)   alpha(n) = [0.1(V+55)] / [1 - exp(-.1(V + 55)]
%       beta(n) = 0.125exp(-.0125(V + 65))
%       alpha(m) = [.1(V + 40)] / (1 - exp(-.1(V + 40))]
%       beta(m) = 4exp(-.0556(V + 65))
%       alpha(h) = .07exp(-.05(V + 65))
%       beta(h) = 1 / [1 + exp(-.1(V + 35))]
%
% In general, equations (2) and (3) can be modeled via the Euler approx.
% or using the ODE solution.
% ------------------------------------------------------------------------
%    
% NOTE: updating of equations (2) and (3) should alternate time steps 
% in order to preserve accuracy to second order in dt.
%
%
%----------------------------------------------------------
% based on code from summer 2003 with Rajesh Rao
% Brian Lundstrom, 6/04
% 10/04, 11/05

% Parameters:
   cm = 10; %nF/mm2
   Gl = 0.003; %mS/mm2
   Gk = 0.36; %mS/mm2
   Gna = 1.2; %mS/mm2
   El = -54.387; %mV
   Ek = -77; %mV
   Ena = 50; %mV

% Initial values:
   Vi = -65; %mV
   mi = 0.0529;
   hi = 0.5961;
   ni = 0.3177;

% TIME
steps = SampleRate; %time steps per ms
t = length(Stim)/steps;
time = t*steps; %ms*steps
dt = 2/steps; %ms, i.e. updating every 2 time steps

% ZERO MEAN NOISE
%randn('state',sum(100*clock));

% EXPONENTIAL FILTER
tau = 2; %time constant (msec)
tmp = tau * SampleRate; % time steps per tau
filt = exp(-(0:(7*tmp))/(tmp));
filt = filt/norm(filt);
ie = conv2(Stim,filt,'full');

% GAUSSIAN FILTER
%half_width = 1.5*SampleRate;
%gvar = 1.3; %variance (msec)
%filt = exp(-(-half_width:half_width).^2/(2*gvar*SampleRate/1000));
%ie = conv(stim,filt);

% Create table for V, n, m, and h
%   where time increases with increasing column number.
% V is updated during odd column numbers.
% n, m, and h are updated during even column numbers.

V = zeros(1,time);
n = zeros(1,time);
m = zeros(1,time);
h = zeros(1,time);

im = zeros(1,time);
V(1) = Vi;
n(1) = ni;
m(1) = mi;
h(1) = hi;

for c = 2:time
    
    if rem(c,2) == 0 %update V during even time steps
   
    % Equations (1) and (2), Euler approx.
    im(c-1) = 1000*(Gl*(V(c-1) - El) + Gk*n(c-1)^4*(V(c-1) - Ek) + ...
        Gna*m(c-1)^3*h(c-1)*(V(c-1) - Ena));
    V(c) = V(c-1) - dt*(im(c-1)/cm - ie(c-1)/cm);
    
    n(c) = n(c-1);
    m(c) = m(c-1);
    h(c) = h(c-1);
    im(c) = im(c-1);
     
    else % update n,m,h during odd time steps
    
    % NOTE: there is a small chance that V = -55.0 or -40.0, in which
    % case the model will break down due to the singularity for 
    % alpha_n and alpha_m, respectively
    alpha_n = (0.01 * (V(c-1) + 55)) / (1 - exp(-.1*(V(c-1) + 55)));
    beta_n = 0.125 * exp(-.0125*(V(c-1) + 65));
    alpha_m = (0.1 * (V(c-1) + 40)) / (1 - exp(-.1*(V(c-1) + 40)));
    beta_m = 4 * exp(-.0556*(V(c-1) + 65));
    alpha_h = .07 * exp(-.05*(V(c-1) + 65));
    beta_h = 1 / (1 + exp(-.1*(V(c-1) + 35)));
    
    % Equation (3), Euler approx.
    n(c) = n(c-1) + dt*( alpha_n*(1-n(c-1)) - beta_n*n(c-1));
    m(c) = m(c-1) + dt*( alpha_m*(1-m(c-1)) - beta_m*m(c-1));
    h(c) = h(c-1) + dt*( alpha_h*(1-h(c-1)) - beta_h*h(c-1));
      
    V(c) = V(c-1);
    
    end
end

y_eul(:,1) = V; y_eul(:,2) = n; y_eul(:,3) = m; y_eul(:,4) = h;
t_eul = (1:time)./steps;

% remove the odd time steps (V changes on even, nmh change on odd)
y = y_eul(2:2:length(y_eul),:);
t = t_eul(2:2:length(t_eul));
ie = ie(2:2:length(ie));

%
figure;
subplot(311)
plot(t,ie(1:length(t)))
ylabel('Current (nA/mm^2)')
title('Hodgkin and Huxley Neuron, approx. via Euler Method')
subplot(312)
plot(t,y(:,1),'b-')
ylabel('V (mV)')
subplot(313);
plot(t,y(:,2),t,y(:,3),t,y(:,4));
ylabel('Gating Variables')
xlabel('msec')
legend('n','m','h')
%}
 