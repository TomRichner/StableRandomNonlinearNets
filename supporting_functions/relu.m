function [y] = relu(x)
% halfwave rectifier

y = x;
y(x<0) = 0;