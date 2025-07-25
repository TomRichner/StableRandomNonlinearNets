function [tout yout] = hhrk4(stim, tspan, dt, y)
% [tout yout] = hhrk4(stim, tspan, dt, y) 
%
% Hodgekin-Huxley model neuron implemented 
% with a 4th Order Runge-Kutta Integrator
%
% STIM = injected current in nA/mm^2
% TSPAN = 2 element vector with begining and end time in msec, e.g. [0 1000];
% DT = time step in msec, e.g. 0.05 msec
% Y = initial values for voltage and gating variables n, m, and h,
%     which are commonly [0 0.3177 0.0529 0.5961].
%
% TOUT = 1xn time vector in msec
% YOUT = 4xn vector of voltage (rest = 0 mV) and n,m, and h gating variables
%
%   NOTE: STIM vector MUST be at least 2x length of your expected output,
%   i.e. TOUT and YOUT
%   This is because the RK4 integrator samples at dt and dt/2.
%
% ------------------------------------------------------------------------
% Brian Lundstrom, Autumn 2004
