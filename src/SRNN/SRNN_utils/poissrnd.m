function R = poissrnd(lambda)
%POISSRND Generate Poisson random numbers without Statistics Toolbox.
%   R = poissrnd(lambda) returns an array the same size as lambda, with
%   Poisson-distributed random numbers generated using Knuth's algorithm.
%
%   Input:
%     lambda - array of non-negative expected counts.
%   Output:
%     R      - array of the same size as lambda, containing Poisson random samples.

% Preallocate output
R = zeros(size(lambda));

% Flatten lambda for iteration
tmp = lambda(:);
R_flat = zeros(size(tmp));

for idx = 1:numel(tmp)
    lam = tmp(idx);
    if lam < 0
        error('poissrnd:NegativeLambda','Lambda must be non-negative.');
    end
    % Knuth's algorithm
    L = exp(-lam);
    k = 0;
    p = 1;
    while p > L
        k = k + 1;
        p = p * rand();
    end
    R_flat(idx) = k - 1;
end

% Reshape back to original size
R = reshape(R_flat, size(lambda));
end 