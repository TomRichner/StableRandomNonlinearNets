function [y] = satu(x,max_x)
% halfwave rectifier

y = x;
y(x>max_x) = max_x;