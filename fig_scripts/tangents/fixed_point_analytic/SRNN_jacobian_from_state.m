function [J, lambda_max, gamma] = SRNN_jacobian_from_state( ...
    M, isE, cSFA_vec, tau_a, tau_b, tau_d, tau_STD, ...
    u, A, b, r, F_STD, Fa)
% Build the continuous-time Jacobian of your SRNN at an arbitrary state
% (u, A, b, r) and return its largest-real-part eigenvalue.
%
% Conventions match your SRNN_NL:
%   r = ReLU(u - cSFA .* sum_k a_k)                (SFA only where cSFA>0)
%   tau_a(k) * a_k' = -a_k + Fa(k) * r             (E only)
%   tau_b * b'      = 1 - b - (b .* r)/tau_STD     (STD only where F_STD>0)
%   tau_d * u'      = -u + DC + M * (r .* b_full), b_full=1 on non-STD
%
% Inputs:
%   M        n×n
%   isE      n×1 logical
%   cSFA_vec n×1 (per-neuron)
%   tau_a    K×1
%   tau_b    scalar
%   tau_d    scalar
%   tau_STD  scalar
%   u        n×1 (current dendrites)
%   A        n×K (SFA states; zeros on non-SFA neurons)
%   b        n×1 (STD states; 1 on non-STD neurons)
%   r        n×1 (current rates)
%   F_STD    n×1 (1 where STD is on, else 0)
%   Fa       K×1 (SFA gains; default = ones(K,1))
%
% Outputs:
%   J          sparse Jacobian
%   lambda_max largest real part of eigenvalues(J)  [1/s]
%   gamma      n×1 ReLU slope mask at this state

n = size(M,1); isE = logical(isE(:));
K = numel(tau_a); tau_a = tau_a(:);
cSFA_vec = cSFA_vec(:);
if nargin < 13 || isempty(Fa), Fa = ones(K,1); else, Fa = Fa(:); end
F_STD = F_STD(:);

% ReLU slope mask at this state (match your activation)
u_eff = u - cSFA_vec .* sum(A,2);
gamma = double(u_eff > 0);                 % d(ReLU)/du
Gam   = spdiags(gamma,0,n,n);

% Indices/selectors
E_idx     = find(isE);
isSTD     = isE & (F_STD > 0);
Estd_idx  = find(isSTD);
nE        = numel(E_idx);
nEstd     = numel(Estd_idx);
P_E       = sparse(E_idx,    1:nE,    1, n, nE);      % n×nE
P_Estd    = sparse(Estd_idx, 1:nEstd, 1, n, nEstd);   % n×nEstd
S = P_Estd.' * P_E;  % maps E -> Estd

% Per-neuron diagonals
bfull = ones(n,1); bfull(isE) = b(isE);
cSFA_D  = spdiags(cSFA_vec,0,n,n);
cSFA_E  = spdiags(cSFA_vec(isE),0,nE,nE);
GamE    = spdiags(gamma(isE),0,nE,nE);
bEstd   = b(Estd_idx);
rEstd   = r(Estd_idx);
GamEstd = spdiags(gamma(Estd_idx),0,nEstd,nEstd);

% Block sizes and indexing
Nu = n; Na = nE*K; Nb = nEstd; N = Nu + Na + Nb;
J  = spalloc(N, N, n + nnz(M) + K*nE + K*K*nE + 3*nE + nnz(M(:,Estd_idx)));

iu = 1:n;
ia = cell(K,1);
for k=1:K, ia{k} = n + (k-1)*nE + (1:nE); end
ib = n + K*nE + (1:nEstd);

%% u-row blocks
% J_uu = (-I + M*diag(b.*gamma)) / tau_d
J(iu, iu) = (-speye(n) + M * spdiags(bfull.*gamma,0,n,n)) / tau_d;

% J_u,a_k = (1/tau_d) * M * diag(b) * d r / d a_k
% dr/da_k = -diag(cSFA .* gamma) on SFA neurons (typically E)
Bua = (1/tau_d) * M * ( spdiags(bfull,0,n,n) * (-cSFA_D) * Gam ) * P_E; % n×nE
for k=1:K
    J(iu, ia{k}) = Bua;
end

% J_u,b = (1/tau_d) * M(:,Estd) * diag(r(Estd))
if nEstd>0
    J(iu, ib) = (1/tau_d) * M(:,Estd_idx) * spdiags(rEstd,0,nEstd,nEstd);
end

%% a_k-row blocks (E only)
for k=1:K
    % a_k,u : (Fa(k)/tau_a(k)) * diag(gamma_E) * P_E'
    J(ia{k}, iu) = (Fa(k)/tau_a(k)) * GamE * P_E.';
    % a_k,a_k' : -I/tau_a(k) + (Fa(k)/tau_a(k)) * dr/da_{k'}
    % with dr/da_{k'} = -diag(cSFA_E .* gamma_E)
    for kp=1:K
        blk = (Fa(k)/tau_a(k)) * ( -cSFA_E ) * GamE;  % coupling via r(a)
        if kp==k
            blk = blk - (1/tau_a(k)) * speye(nE);
        end
        J(ia{k}, ia{kp}) = blk;
    end
end

%% b-row blocks (E with STD only)
if nEstd>0
    % b,u : -(1/tau_STD) * diag(b .* gamma) * P_Estd'
    J(ib, iu) = -(1/tau_STD) * GamEstd * spdiags(bEstd,0,nEstd,nEstd) * P_Estd.';
    % b,a_k : + (1/tau_STD) * diag(b .* gamma .* cSFA)  (restricted to Estd)
    C = (1/tau_STD) * GamEstd * spdiags(bEstd,0,nEstd,nEstd) * (S * cSFA_E);
    for k=1:K
        J(ib, ia{k}) = C;
    end
    % b,b : -diag(1/tau_b + r/tau_STD)
    J(ib, ib) = -spdiags( (1/tau_b) + (rEstd/tau_STD), 0, nEstd, nEstd );
end

%% Largest-real-part eigenvalue
try
    lambda_max = max(real(eigs(J,1,'lr')));
catch
    lambda_max = max(real(eig(full(J))));
end
end
