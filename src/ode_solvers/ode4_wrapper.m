function [t_ode, Y] = ode4_wrapper(odefun, tspan, y0, varargin)
%ODE4_WRAPPER Wraps ode4 to provide an output compatible with ode45.
%   [T,Y] = ODE4_WRAPPER(ODEFUN,TSPAN,Y0,...) calls the non-adaptive 
%   ODE solver ODE4 and returns the time vector TSPAN (as a column vector T)
%   and the solution array Y. This makes its syntax similar to MATLAB's
%   adaptive solvers like ODE45.
%
%   The underlying ODE4 function must be in the MATLAB path.
%
%   Example:
%       tspan = 0:0.1:20;
%       [t,y] = ode4_wrapper(@vdp1,tspan,[2 0]);  
%       plot(t,y(:,1));
%
%   See also ODE4, ODE45.

% Call the original ode4 function
Y_internal = ode4(odefun, tspan, y0);

% Prepare outputs similar to ode45
% tspan is used as the time vector, ensure it's a column vector
t_ode = tspan(:); 
Y = Y_internal;

end