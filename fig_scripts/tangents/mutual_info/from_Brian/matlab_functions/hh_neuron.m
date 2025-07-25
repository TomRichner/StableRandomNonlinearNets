% Hodgkin and Huxley model of neuronal membrane
%
% Brian Lundstrom, Autumn 2004

clear

%__________________________________________________________
%
% parameters
stim_mean = 50; % nA/mm^2
stim_var = 100;
tspan = [0 3000]; % msec
dt = 0.05; % msec
%
%__________________________________________________________

steps = dt/2; % use dt/2 as step size since this is required by RK4 solver
time = tspan(1):steps:tspan(2);

% make stimulus; STIMBIT2 is used when concatenating stimulus vectors
% NOTE: stimulus length needs to be sampled at 2x of expected output since
% RK4 solvers sample at the half time step
stim_length = 2*tspan(2)/1000;
[stim stimbit2] = make_stim(stim_length,stim_mean,stim_var,20000,[]);

% use Hodgkin-Huxley Runge-Kutta 4th order solver with appropriate initial
% values
[t y] = hhrk4(stim,tspan,dt,[0 0.3177 0.0529 0.5961]);

% plot results
clf
x = time;
subplot(3,1,1)
plot(time,stim(1:length(time)))
title('Hodgkin Huxley Neuron, approx. via RK4')
ylabel('Stimulus')
subplot(3,1,2);
plot(t,y(1,:))
ylabel('Membrane Potential (mV)')
subplot(3,1,3);
plot(t,y(2,:),t,y(3,:),t,y(4,:));
xlabel('msec')
ylabel('Gating Variables')
legend('n','m','h');