function [J, lambda_max, gamma] = SRNN_jacobian_from_state( ...
    M, isE, cSFA_vec, tau_a, tau_b, tau_d, tau_STD, ...
    u, A, b, r, F_STD)
% Build continuous-time Jacobian at the provided state (u,A,b,r),
% using your SRNN_NL equations and ReLU derivative mask.
% Works at a fixed point OR any arbitrary state (e.g., final state of a sim).
%
% A is n×K with zeros on non-SFA neurons; b is n×1 with ones where STD off.

n   = size(M,1); isE = logical(isE(:)); isI = ~isE;
K   = numel(tau_a); tau_a = tau_a(:);
cSFA_vec = cSFA_vec(:);
F_STD = F_STD(:);
nE  = nnz(isE);

% ReLU slope mask at this state
u_eff = u - cSFA_vec .* sum(A,2);   % same as SRNN_NL
gamma = double(u_eff > 0);          % d(ReLU)/du

% Indexing helpers
E_idx    = find(isE);
Estd_idx = find(isE & (F_STD>0));    % STD-active E-neurons
nEstd    = numel(Estd_idx);

% Selectors
P_E    = sparse(E_idx,    1:nE,    1, n, nE);     % n×nE
P_Estd = sparse(Estd_idx, 1:nEstd, 1, n, nEstd);  % n×nEstd
S = P_Estd.' * P_E;  % (nEstd×nE) maps from E to E_std subset

% Full b vector (1 where STD off)
bfull = ones(n,1); bfull(isE) = b(isE);

% Block sizes
Nu   = n;
Na   = nE*K;
Nb   = nEstd;
Ntot = Nu + Na + Nb;

J = spalloc(Ntot, Ntot,  ...
     n + nnz(M) + K*nE + K*K*nE + 3*nE + nnz(M(:,Estd_idx)));

idx_u = 1:n;
idx_a = cell(K,1);
for k=1:K, idx_a{k} = n + (k-1)*nE + (1:nE); end
idx_b = n + K*nE + (1:nEstd);

% Shorthands
Gam    = spdiags(gamma,0,n,n);
GamE   = spdiags(gamma(isE),0,nE,nE);
cSFA_D = spdiags(cSFA_vec,0,n,n);
cSFA_E = spdiags(cSFA_vec(isE),0,nE,nE);

%% J_uu = (-I + M*diag(b.*gamma)) / tau_d
J(idx_u, idx_u) = (-speye(n) + M * spdiags(bfull.*gamma,0,n,n)) / tau_d;

%% J_u,a_k = (1/tau_d) * M * diag(b) * d r/ d a_k
% dr/da_k = -diag(cSFA .* gamma) on SFA neurons (typically E)
Bua = (1/tau_d) * M * ( spdiags(bfull,0,n,n) * (-cSFA_D) * Gam ) * P_E; % n×nE
for k=1:K
    J(idx_u, idx_a{k}) = Bua;
end

%% J_u,b = (1/tau_d)*M(:,Estd) * diag(r(Estd))
if nEstd>0
    J(idx_u, idx_b) = (1/tau_d) * M(:,Estd_idx) * spdiags(r(Estd_idx),0,nEstd,nEstd);
end

%% SFA rows
for k=1:K
    % a_k,u : (1/tau_a(k)) * diag(gamma_E) * P_E'
    J(idx_a{k}, idx_u) = (1/tau_a(k)) * GamE * P_E.';
    % a_k,a_k' : -I/tau_a(k) + (1/tau_a(k)) * dr/da_{k'}
    % dr/da_{k'} = -diag(cSFA_E * gamma_E)
    for kp=1:K
        blk = (1/tau_a(k)) * (-cSFA_E) * GamE;  % coupling from sum(a)
        if kp==k, blk = blk - (1/tau_a(k))*speye(nE); end
        J(idx_a{k}, idx_a{kp}) = blk;
    end
end

%% STD rows (only on Estd_idx)
if nEstd>0
    % b,u : -(1/tau_STD) * diag(b.*gamma) * P_Estd'
    J(idx_b, idx_u) = -(1/tau_STD) * spdiags(b(Estd_idx).*gamma(Estd_idx),0,nEstd,nEstd) * P_Estd.';
    % b,a_k : +(1/tau_STD) * diag(b.*gamma .* cSFA) * S   (S selects the Estd subset)
    C = (1/tau_STD) * spdiags(b(Estd_idx).*gamma(Estd_idx),0,nEstd,nEstd) * (S * cSFA_E);
    for k=1:K
        J(idx_b, idx_a{k}) = C;
    end
    % b,b : -diag(1/tau_b + r/τ_STD)
    J(idx_b, idx_b) = -spdiags( (1/tau_b) + (r(Estd_idx)/tau_STD), 0, nEstd, nEstd );
end

%% Spectral abscissa
try
    lambda_max = max(real(eigs(J,1,'lr')));
catch
    lambda_max = max(real(eig(full(J))));
end
end
