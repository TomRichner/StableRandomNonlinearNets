function [u, A, b, r, info] = SRNN_fixed_point( ...
    M, isE, cSFA_vec, tau_a, tau_b, tau_d, tau_STD, DC, F_STD, opts)
% Solve the algebraic steady state of your SRNN for a constant DC input.
% Matches SRNN_NL conventions exactly:
%   - SFA only where cSFA_vec(i)~=0 (typically E only)
%   - STD only where F_STD(i)>0 (typically E only)
%
% Inputs:
%   M          n×n
%   isE        n×1 logical
%   cSFA_vec   n×1 vector (per-neuron SFA coefficient; e.g., 0.5/K on E, 0 on I)
%   tau_a      K×1 vector of SFA time constants (E only)
%   tau_b      scalar STD recovery time constant
%   tau_d      scalar dendrite time constant
%   tau_STD    scalar depletion time constant (in db/dt term)
%   DC         scalar or n×1 constant input
%   F_STD      n×1 mask (1 where STD is active, 0 else). Pass [] to default to E-only.
%   opts.tol, opts.maxit, opts.alpha, opts.Fa (K×1 gains, default 1)
%
% Outputs:
%   u  (n×1), A (n×K; zeros off SFA neurons), b (n×1; ones where STD off), r (n×1)
%   info struct

n = size(M,1); isE = logical(isE(:)); isI = ~isE;
K = numel(tau_a); tau_a = tau_a(:);
cSFA_vec = cSFA_vec(:);
if isempty(F_STD), F_STD = double(isE); else, F_STD = F_STD(:); end

% Options
if nargin < 10 || isempty(opts), opts = struct; end
tol     = get_opt(opts,'tol',1e-11);
maxit   = get_opt(opts,'maxit',5e4);
alpha   = get_opt(opts,'alpha',0.6);
Fa      = get_opt(opts,'Fa',ones(K,1));
verbose = get_opt(opts,'verbose',true);

% DC vector
if isscalar(DC), DC = DC*ones(n,1); else, DC = DC(:); end

% Per-neuron total SFA gain Phi_i = cSFA_i * sum_k Fa_k  (typically = 0.5 on E)
Phi = cSFA_vec * sum(Fa);

% Init
u = DC; r = zeros(n,1);
b = ones(n,1);
A = zeros(n,K);

for it = 1:maxit
    % STD steady-state b from current r, only where STD is ON
    r_pos = max(r,0);
    b_new = ones(n,1);
    idx_std = F_STD>0;
    b_new(idx_std) = tau_STD ./ (tau_STD + tau_b * r_pos(idx_std));  % in [0,1]

    % Dendrites at steady state: u = DC + M*(r .* b_eff)
    s = r_pos .* b_new;
    u_new = DC + M * s;

    % Rates at steady state given u_new (piecewise linear)
    r_new = zeros(n,1);
    % On neurons with SFA (typically E): r = ReLU(u/(1+Phi_i))
    idx_sfa = (Phi > 0);
    r_new(idx_sfa) = max(0, u_new(idx_sfa) ./ (1 + Phi(idx_sfa)));
    % Others (typically I): pure ReLU
    idx_no_sfa = ~idx_sfa;
    r_new(idx_no_sfa) = max(0, u_new(idx_no_sfa));

    % Damped update
    u = (1-alpha)*u + alpha*u_new;
    r = (1-alpha)*r + alpha*r_new;

    % Convergence (∞-norm on both)
    res = max(norm(u_new - u, inf), norm(r_new - r, inf));
    if res < tol, break; end
end

% Final b and A
b = ones(n,1);
b(idx_std) = tau_STD ./ (tau_STD + tau_b * max(r(idx_std),0));
A = zeros(n,K);
if K>0
    A(idx_sfa,:) = (max(r(idx_sfa),0)) * (Fa(:).');  % A_i,k = Fa_k * r_i
end

info = struct('converged', it<maxit, 'niters', it, 'residual', res, ...
              'tol', tol, 'alpha', alpha);
if verbose
    fprintf('SRNN_fixed_point: %s in %d iters (residual %.3g)\n', ...
        tern(info.converged,'converged','stopped'), it, res);
end
end

% helpers
function v = get_opt(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function t = tern(c,a,b), if c, t=a; else, t=b; end, end
