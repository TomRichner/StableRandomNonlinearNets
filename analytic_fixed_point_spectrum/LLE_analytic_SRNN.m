%% -----------------------------------------------------------
%  File:   LLE_analytic_SRNN.m
%
%  Analytic fixed‑point and Lyapunov calculations for the
%  SRNN rate model with ReLU, multi‑τ SFA and single‑τ STD.
%
%  Assumptions:
%     • identical SFA / STD parameters for all excitatory neurons
%     • n_b = 1  (single STD time‑constant)
%     • all neurons receive the same constant DC input u_ex
%     • all excitatory neurons are active at the fixed point
%       (r0_E > 0) ; inhibitory neurons may be silent or active
%
%  -----------------------------------------------------------

clear;  clc;

%% ------------------ 1.   USER‑DEFINED PARAMETERS -----------------------

n     = 10;               % total neurons
fracE = 0.7;               % E fraction; remainder are I
nE    = round(fracE*n);
nI    = n - nE;

% Connectivity ------------------------------------------------------------
g     = 0.5;               % weight scale (σ of Gaussian entries)
density = 0.55;             % connection probability
% rng(42);
M = g/sqrt(n) * (randn(n) .* (rand(n)<density));   % dense + sparse mask
% enforce Dale's law
M(:,1:nE)  =  abs(M(:,1:nE));   % E columns positive
M(:,nE+1:end) = -abs(M(:,nE+1:end));
M(1:n+1:end) = 0;               % zero self-connections

% External drive (DC) -----------------------------------------------------
u_ex = 5 * ones(n,1);   % identical DC to all neurons

% SFA (excitatory only) ---------------------------------------------------
n_a   = 3;
tau_a = logspace(log10(0.3),log10(15),n_a)';  % column!
c_SFA = ones(n,1);          % strength; zeroed for I below
c_SFA(nE+1:end) = 0;        % no SFA on I cells

% STD (single τ on E only) -------------------------------------------------
tau_b  = 2;                 % s
tau_STD = 0.5;              % s
F_STD  = ones(n,1);
F_STD(nE+1:end) = 0;        % no STD on I cells

tau_d  = 0.025;             % dendritic low‑pass

tol_fp = 1e-9;             % fixed‑point tolerance
maxIt  = 10000;

%% ------------------ 2.   FIXED‑POINT ITERATION -------------------------

% helper: gamma = F_STD * tau_b / tau_STD (only non‑zero on E)
gamma = F_STD .* (tau_b/tau_STD);

alpha  = 0.1;          % 0<alpha≤1  (smaller ⇒ safer, slower)
r      = max(0,u_ex);  % initial guess
for it = 1:maxIt
    p = r ./ (1 + gamma .* r);
    u_d = u_ex + M*p;
    u_eff = u_d - c_SFA .* n_a .* r;
    r_raw = max(0, u_e
    ff);
    r_new = (1-alpha)*r + alpha*r_raw;   % Damped update
    
    if norm(r_new - r,inf) < tol_fp
        r = r_new;  break
    end
    
    if it==maxIt,  error('Fixed point did not converge'); end

    r = r_new;
end
fprintf('Fixed‑point solved in %d iterations\n',it);

r0 = r;
p0 = r0 ./ (1 + gamma .* r0);
h0 = 1 ./ (1 + gamma .* r0);
b0 = h0;                     % single‑τ => b0 == h0

% sanity: print max firing rate
fprintf('max(r0) = %.4g Hz\n',max(r0));

%% ------------------ 3.   BUILD POLYNOMIAL COEFFICIENTS -----------------

% Convenience: same parameters for all *active* E neurons
idxActive = find(r0>0);      % set 𝔄
if any(c_SFA(idxActive)~=c_SFA(idxActive(1))) || ...
   any(gamma(idxActive)   ~=gamma(idxActive(1)))
    warning('Heterogeneous adaptation parameters – analytic block‑diagonal derivation is needed.');
end
cS  = c_SFA(idxActive(1));
rS  = r0(idxActive(1));
bS  = b0(idxActive(1));
hS  = h0(idxActive(1));
gammaS = gamma(idxActive(1));

% ---- A(λ)   (SFA) -------------------------------------------------------
DenA = 1;                       % coefficients, highest power first
for k = 1:n_a
    DenA = conv(DenA,[tau_a(k) 1]);     % multiply (τ_a λ + 1)
end

NumerA = DenA;                  % start with "1"
for k = 1:n_a
    % divide by (τ_a λ + 1) to get prod_{j≠k}
    Partial = deconv(DenA,[tau_a(k) 1]);
    NumerA = polyadd(NumerA, cS*Partial);   %#ok<*NOPRT>
end

% ---- B(λ)   (STD) -------------------------------------------------------
K   = F_STD(idxActive(1))*rS/(bS*tau_STD);   % dimensionless
DenB   = [1 1/tau_b];              % λ + 1/τ_b
NumerB = [1 1/tau_b+K];            % λ + 1/τ_b + K

% ---- C(λ) = (λ τ_d + 1) -----------------------------------------------
Cpoly = [tau_d 1];

% helper for poly multiplication
ABC_num = conv(conv(Cpoly,NumerA),NumerB);   % (λτ_d+1)·NumerA·NumerB
ABD_den = conv(DenA,DenB);                  % DenA·DenB   (degree n_a+1)

%% ------------------ 4.   EIGEN‑LOOP  &  LLE ----------------------------

eigM = eig(M);           % μ_k
LLE  = -Inf;

for k = 1:n
    mu = eigM(k);
    % polynomial: ABC_num  -  mu*hS*ABD_den  = 0
    RHS = mu*hS*ABD_den;
    % pad shorter vector with zeros to equal length
    len = max(numel(ABC_num),numel(RHS));
    p_left  = [zeros(1,len-numel(ABC_num))  ABC_num];
    p_right = [zeros(1,len-numel(RHS))      RHS];
    polyEq  = p_left - p_right;          % coefficients of λ^N ... λ^0
    
    lambdaRoots = roots(polyEq);
    LLE = max(LLE, max(real(lambdaRoots)));
end
fprintf('Analytic largest Lyapunov exponent  Λ_max = %+8.5f  1/s\n',LLE);

%% ------------------ 5.   HELPER (polyadd) ------------------------------
function s = polyadd(a,b)
% add two polynomial coefficient vectors (highest power first)
if numel(a)<numel(b)
    a = [zeros(1,numel(b)-numel(a)) a];
elseif numel(b)<numel(a)
    b = [zeros(1,numel(a)-numel(b)) b];
end
s = a+b;
end
