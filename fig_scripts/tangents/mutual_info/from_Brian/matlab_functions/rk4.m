function [T,Y] = rk4(func_ty, tspan, dt, y)
%[T,Y] = rk4(func_ty, tspan, dt, y)
% Fourth-order Runge-Kutta
%
% call with: 
% [T, Y] = rk4(@hheval,[5 5 5],[0;0;0;0]);
%
% Brian Lundstrom, Autumn 2004, 1/08

dt2 = dt/2.0;
tt = tspan(1):dt:tspan(2)-dt;
T = zeros(1,length(tt));
Y = zeros(length(y),length(tt));
T(1) = 0;
Y(:,1) = y;

for i = 1:length(tt)
    t = tt(i);
        
    k1 = dt*feval(func_ty, t, y);
    k2 = dt*feval(func_ty, t+dt2, y+0.5.*k1);
    k3 = dt*feval(func_ty, t+dt2, y+0.5.*k2);
    k4 = dt*feval(func_ty, t+dt, y+k3);
    y = y + 1/6.*(k1 + k4) + 1/3.*(k2 + k3);
    
    T(i+1) = t+dt;
    Y(:,i+1) = y;
    
end 
